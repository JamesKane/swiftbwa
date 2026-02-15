#if canImport(Metal)
import Metal
import BWACore

/// A single banded global SW + traceback task to dispatch to the GPU.
public struct GlobalSWTask: Sendable {
    public let query: [UInt8]      // 2-bit encoded query segment
    public let target: [UInt8]     // 2-bit encoded ref segment
    public let trueScore: Int32    // expected score for bandwidth retry
    public let initialW: Int32     // initial bandwidth
    public let scoring: ScoringParameters

    public init(query: [UInt8], target: [UInt8], trueScore: Int32,
                initialW: Int32, scoring: ScoringParameters) {
        self.query = query
        self.target = target
        self.trueScore = trueScore
        self.initialW = initialW
        self.scoring = scoring
    }
}

/// Result from GPU global SW: score + packed CIGAR ops.
public struct GlobalSWResult: Sendable {
    public var score: Int32
    public var cigar: [UInt32]     // packed CIGAR ops (length << 4 | op)

    public init(score: Int32, cigar: [UInt32]) {
        self.score = score
        self.cigar = cigar
    }
}

private let MAX_CIGAR_OPS = 512

/// Batch dispatcher for banded global SW + CIGAR on Metal GPU.
public struct GlobalSWDispatcher: Sendable {

    /// Dispatch a batch of global SW tasks to the GPU synchronously.
    /// Returns one GlobalSWResult per task.
    public static func dispatchBatch(
        tasks: [GlobalSWTask],
        engine: MetalSWEngine
    ) -> [GlobalSWResult] {
        guard !tasks.isEmpty, let pipeline = engine.globalSWPipeline else {
            return []
        }

        let taskCount = tasks.count

        // Pack queries into contiguous buffer with offsets
        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []

        var allTargets: [UInt8] = []
        var targetOffsets: [UInt32] = []
        var targetLengths: [UInt16] = []

        var trueScoreValues: [Int16] = []
        var initialWValues: [Int16] = []

        var wsOffsets: [UInt32] = []
        var totalWorkspaceBytes: UInt32 = 0

        for task in tasks {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(task.query.count))
            allQueries.append(contentsOf: task.query)

            targetOffsets.append(UInt32(allTargets.count))
            targetLengths.append(UInt16(task.target.count))
            allTargets.append(contentsOf: task.target)

            trueScoreValues.append(Int16(clamping: task.trueScore))
            initialWValues.append(Int16(clamping: task.initialW))

            let qlen = task.query.count
            let tlen = task.target.count

            // Workspace: qp (5*qlen*4) + H ((qlen+1)*4) + E ((qlen+1)*4) + z (tlen*nCol)
            // For z, we need worst-case nCol after retries (w can grow up to 8x).
            // Max w after 3 retries = initialW * 8, but capped by qlen.
            // nCol = min(qlen, 2*w+1) so worst case nCol = qlen.
            let maxNCol = qlen  // worst case: full bandwidth
            let wsBytes = UInt32(5 * qlen * 4 + (qlen + 1) * 4 + (qlen + 1) * 4 + tlen * maxNCol)

            // Align to 16 bytes
            let aligned = (totalWorkspaceBytes + 15) & ~15
            wsOffsets.append(aligned)
            totalWorkspaceBytes = aligned + wsBytes
        }

        // Use first task's scoring params (all tasks in a batch share scoring)
        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let params: [Int16] = [
            Int16(scoring.gapOpenPenalty),
            Int16(scoring.gapExtendPenalty),
            Int16(scoring.gapOpenPenaltyDeletion),
            Int16(scoring.gapExtendPenaltyDeletion)
        ]

        let pool = engine.bufferPool

        // Create Metal buffers
        guard let queriesBuf = pool.acquire(minSize: max(allQueries.count, 1)),
              let queryOffsetsBuf = pool.acquire(minSize: queryOffsets.count * 4),
              let queryLengthsBuf = pool.acquire(minSize: queryLengths.count * 2),
              let targetsBuf = pool.acquire(minSize: max(allTargets.count, 1)),
              let targetOffsetsBuf = pool.acquire(minSize: targetOffsets.count * 4),
              let targetLengthsBuf = pool.acquire(minSize: targetLengths.count * 2),
              let matBuf = pool.acquire(minSize: mat.count),
              let paramsBuf = pool.acquire(minSize: params.count * 2),
              let trueScoresBuf = pool.acquire(minSize: trueScoreValues.count * 2),
              let initialWBuf = pool.acquire(minSize: initialWValues.count * 2),
              let cigarOutBuf = pool.acquire(minSize: taskCount * MAX_CIGAR_OPS * 4),
              let cigarLengthsBuf = pool.acquire(minSize: taskCount * 4),
              let scoresBuf = pool.acquire(minSize: taskCount * 4),
              let wsBuf = pool.acquire(minSize: max(Int(totalWorkspaceBytes), 1)),
              let wsOffsetsBuf = pool.acquire(minSize: wsOffsets.count * 4)
        else {
            return []
        }

        let allBuffers = [queriesBuf, queryOffsetsBuf, queryLengthsBuf,
                          targetsBuf, targetOffsetsBuf, targetLengthsBuf,
                          matBuf, paramsBuf, trueScoresBuf, initialWBuf,
                          cigarOutBuf, cigarLengthsBuf, scoresBuf, wsBuf, wsOffsetsBuf]

        // Copy data into buffers
        if !allQueries.isEmpty {
            memcpy(queriesBuf.contents(), allQueries, allQueries.count)
        }
        memcpy(queryOffsetsBuf.contents(), queryOffsets, queryOffsets.count * 4)
        memcpy(queryLengthsBuf.contents(), queryLengths, queryLengths.count * 2)
        if !allTargets.isEmpty {
            memcpy(targetsBuf.contents(), allTargets, allTargets.count)
        }
        memcpy(targetOffsetsBuf.contents(), targetOffsets, targetOffsets.count * 4)
        memcpy(targetLengthsBuf.contents(), targetLengths, targetLengths.count * 2)
        memcpy(matBuf.contents(), mat, mat.count)
        memcpy(paramsBuf.contents(), params, params.count * 2)
        memcpy(trueScoresBuf.contents(), trueScoreValues, trueScoreValues.count * 2)
        memcpy(initialWBuf.contents(), initialWValues, initialWValues.count * 2)
        memcpy(wsOffsetsBuf.contents(), wsOffsets, wsOffsets.count * 4)

        // Encode and dispatch
        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(allBuffers)
            return []
        }

        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(queriesBuf, offset: 0, index: 0)
        encoder.setBuffer(queryOffsetsBuf, offset: 0, index: 1)
        encoder.setBuffer(queryLengthsBuf, offset: 0, index: 2)
        encoder.setBuffer(targetsBuf, offset: 0, index: 3)
        encoder.setBuffer(targetOffsetsBuf, offset: 0, index: 4)
        encoder.setBuffer(targetLengthsBuf, offset: 0, index: 5)
        encoder.setBuffer(matBuf, offset: 0, index: 6)
        encoder.setBuffer(paramsBuf, offset: 0, index: 7)
        encoder.setBuffer(trueScoresBuf, offset: 0, index: 8)
        encoder.setBuffer(initialWBuf, offset: 0, index: 9)
        encoder.setBuffer(cigarOutBuf, offset: 0, index: 10)
        encoder.setBuffer(cigarLengthsBuf, offset: 0, index: 11)
        encoder.setBuffer(scoresBuf, offset: 0, index: 12)
        encoder.setBuffer(wsBuf, offset: 0, index: 13)
        encoder.setBuffer(wsOffsetsBuf, offset: 0, index: 14)

        let threadgroupSize = min(64, pipeline.maxTotalThreadsPerThreadgroup)
        let gridSize = MTLSize(width: taskCount, height: 1, depth: 1)
        let tgSize = MTLSize(width: threadgroupSize, height: 1, depth: 1)
        encoder.dispatchThreads(gridSize, threadsPerThreadgroup: tgSize)
        encoder.endEncoding()

        cmdBuf.commit()
        cmdBuf.waitUntilCompleted()

        // Read results
        var results: [GlobalSWResult] = []
        results.reserveCapacity(taskCount)
        let cigarPtr = cigarOutBuf.contents().bindMemory(to: Int32.self, capacity: taskCount * MAX_CIGAR_OPS)
        let cigarLenPtr = cigarLengthsBuf.contents().bindMemory(to: Int32.self, capacity: taskCount)
        let scorePtr = scoresBuf.contents().bindMemory(to: Int32.self, capacity: taskCount)

        for i in 0..<taskCount {
            let nOps = Int(cigarLenPtr[i])
            let score = scorePtr[i]
            var cigar: [UInt32] = []
            cigar.reserveCapacity(nOps)
            let base = i * MAX_CIGAR_OPS
            for j in 0..<nOps {
                cigar.append(UInt32(bitPattern: cigarPtr[base + j]))
            }
            results.append(GlobalSWResult(score: score, cigar: cigar))
        }

        pool.release(allBuffers)
        return results
    }
}
#endif
