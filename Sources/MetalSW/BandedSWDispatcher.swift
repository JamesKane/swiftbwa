#if canImport(Metal)
import Metal
import BWACore

/// A single banded SW task to dispatch to the GPU.
public struct BandedSWTask: Sendable {
    public let query: [UInt8]
    public let target: [UInt8]
    public let h0: Int32
    public let w: Int32
    public let scoring: ScoringParameters

    public init(query: [UInt8], target: [UInt8], h0: Int32, w: Int32, scoring: ScoringParameters) {
        self.query = query
        self.target = target
        self.h0 = h0
        self.w = w
        self.scoring = scoring
    }
}

/// Batch dispatcher for banded Smith-Waterman on Metal GPU.
public struct BandedSWDispatcher: Sendable {

    /// Dispatch a batch of 8-bit banded SW tasks to the GPU.
    /// Returns SWResult for each task. Tasks that overflowed 8-bit have nil result.
    public static func dispatchBatch8(
        tasks: [BandedSWTask],
        engine: MetalSWEngine
    ) -> [SWResult?] {
        guard !tasks.isEmpty else { return [] }

        let taskCount = tasks.count

        // Pack queries into contiguous buffer with offsets
        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []

        var allTargets: [UInt8] = []
        var targetOffsets: [UInt32] = []
        var targetLengths: [UInt16] = []

        var h0Values: [Int16] = []

        // Calculate workspace offsets — each task needs H[qlen] + E[qlen] + profile[5*qlen] bytes
        var wsOffsets: [UInt32] = []
        var totalWorkspaceBytes: UInt32 = 0

        for task in tasks {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(task.query.count))
            allQueries.append(contentsOf: task.query)

            targetOffsets.append(UInt32(allTargets.count))
            targetLengths.append(UInt16(task.target.count))
            allTargets.append(contentsOf: task.target)

            h0Values.append(Int16(clamping: task.h0))

            let qlen = task.query.count
            let wsBytes = UInt32(qlen + qlen + 5 * qlen)  // H + E + profile

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
            Int16(scoring.gapExtendPenaltyDeletion),
            Int16(scoring.zDrop),
            Int16(scoring.mismatchPenalty),
            Int16(tasks[0].w)
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
              let h0Buf = pool.acquire(minSize: h0Values.count * 2),
              let resultsBuf = pool.acquire(minSize: taskCount * 6 * 4),
              let wsBuf = pool.acquire(minSize: max(Int(totalWorkspaceBytes), 1)),
              let wsOffsetsBuf = pool.acquire(minSize: wsOffsets.count * 4)
        else {
            return Array(repeating: nil, count: taskCount)
        }

        let allBuffers = [queriesBuf, queryOffsetsBuf, queryLengthsBuf,
                          targetsBuf, targetOffsetsBuf, targetLengthsBuf,
                          matBuf, paramsBuf, h0Buf, resultsBuf, wsBuf, wsOffsetsBuf]

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
        memcpy(h0Buf.contents(), h0Values, h0Values.count * 2)
        memcpy(wsOffsetsBuf.contents(), wsOffsets, wsOffsets.count * 4)

        // Encode and dispatch
        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(allBuffers)
            return Array(repeating: nil, count: taskCount)
        }

        encoder.setComputePipelineState(engine.bandedSW8Pipeline)
        encoder.setBuffer(queriesBuf, offset: 0, index: 0)
        encoder.setBuffer(queryOffsetsBuf, offset: 0, index: 1)
        encoder.setBuffer(queryLengthsBuf, offset: 0, index: 2)
        encoder.setBuffer(targetsBuf, offset: 0, index: 3)
        encoder.setBuffer(targetOffsetsBuf, offset: 0, index: 4)
        encoder.setBuffer(targetLengthsBuf, offset: 0, index: 5)
        encoder.setBuffer(matBuf, offset: 0, index: 6)
        encoder.setBuffer(paramsBuf, offset: 0, index: 7)
        encoder.setBuffer(h0Buf, offset: 0, index: 8)
        encoder.setBuffer(resultsBuf, offset: 0, index: 9)
        encoder.setBuffer(wsBuf, offset: 0, index: 10)
        encoder.setBuffer(wsOffsetsBuf, offset: 0, index: 11)

        let threadgroupSize = min(64, engine.bandedSW8Pipeline.maxTotalThreadsPerThreadgroup)
        let gridSize = MTLSize(width: taskCount, height: 1, depth: 1)
        let tgSize = MTLSize(width: threadgroupSize, height: 1, depth: 1)
        encoder.dispatchThreads(gridSize, threadsPerThreadgroup: tgSize)
        encoder.endEncoding()

        cmdBuf.commit()
        cmdBuf.waitUntilCompleted()

        // Read results
        var results: [SWResult?] = []
        results.reserveCapacity(taskCount)
        let resultPtr = resultsBuf.contents().bindMemory(to: Int32.self, capacity: taskCount * 6)

        for i in 0..<taskCount {
            let base = i * 6
            let score = resultPtr[base]
            if score <= -2147483647 {
                results.append(nil)  // overflow
            } else {
                results.append(SWResult(
                    score: score,
                    queryEnd: resultPtr[base + 1],
                    targetEnd: resultPtr[base + 2],
                    globalTargetEnd: resultPtr[base + 3],
                    globalScore: resultPtr[base + 4]
                ))
            }
        }

        pool.release(allBuffers)
        return results
    }

    /// Dispatch a batch of 16-bit banded SW tasks to the GPU.
    /// For overflow fallback from 8-bit kernel.
    public static func dispatchBatch16(
        tasks: [BandedSWTask],
        engine: MetalSWEngine
    ) -> [SWResult] {
        guard !tasks.isEmpty, let pipeline = engine.bandedSW16Pipeline else {
            // CPU fallback handled by caller
            return []
        }

        let taskCount = tasks.count

        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []
        var allTargets: [UInt8] = []
        var targetOffsets: [UInt32] = []
        var targetLengths: [UInt16] = []
        var h0Values: [Int16] = []
        var wsOffsets: [UInt32] = []
        var totalWorkspaceBytes: UInt32 = 0

        for task in tasks {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(task.query.count))
            allQueries.append(contentsOf: task.query)

            targetOffsets.append(UInt32(allTargets.count))
            targetLengths.append(UInt16(task.target.count))
            allTargets.append(contentsOf: task.target)

            h0Values.append(Int16(clamping: task.h0))

            // 16-bit: H[qlen] + E[qlen] + profile[5*qlen], each element 2 bytes (short)
            let qlen = task.query.count
            let wsBytes = UInt32((qlen + qlen + 5 * qlen) * 2)

            let aligned = (totalWorkspaceBytes + 15) & ~15
            wsOffsets.append(aligned)
            totalWorkspaceBytes = aligned + wsBytes
        }

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let params: [Int16] = [
            Int16(scoring.gapOpenPenalty),
            Int16(scoring.gapExtendPenalty),
            Int16(scoring.gapOpenPenaltyDeletion),
            Int16(scoring.gapExtendPenaltyDeletion),
            Int16(scoring.zDrop),
            0,  // no bias for signed 16-bit
            Int16(tasks[0].w)
        ]

        let pool = engine.bufferPool

        guard let queriesBuf = pool.acquire(minSize: max(allQueries.count, 1)),
              let queryOffsetsBuf = pool.acquire(minSize: queryOffsets.count * 4),
              let queryLengthsBuf = pool.acquire(minSize: queryLengths.count * 2),
              let targetsBuf = pool.acquire(minSize: max(allTargets.count, 1)),
              let targetOffsetsBuf = pool.acquire(minSize: targetOffsets.count * 4),
              let targetLengthsBuf = pool.acquire(minSize: targetLengths.count * 2),
              let matBuf = pool.acquire(minSize: mat.count),
              let paramsBuf = pool.acquire(minSize: params.count * 2),
              let h0Buf = pool.acquire(minSize: h0Values.count * 2),
              let resultsBuf = pool.acquire(minSize: taskCount * 6 * 4),
              let wsBuf = pool.acquire(minSize: max(Int(totalWorkspaceBytes), 1)),
              let wsOffsetsBuf = pool.acquire(minSize: wsOffsets.count * 4)
        else {
            return []
        }

        let allBuffers = [queriesBuf, queryOffsetsBuf, queryLengthsBuf,
                          targetsBuf, targetOffsetsBuf, targetLengthsBuf,
                          matBuf, paramsBuf, h0Buf, resultsBuf, wsBuf, wsOffsetsBuf]

        if !allQueries.isEmpty { memcpy(queriesBuf.contents(), allQueries, allQueries.count) }
        memcpy(queryOffsetsBuf.contents(), queryOffsets, queryOffsets.count * 4)
        memcpy(queryLengthsBuf.contents(), queryLengths, queryLengths.count * 2)
        if !allTargets.isEmpty { memcpy(targetsBuf.contents(), allTargets, allTargets.count) }
        memcpy(targetOffsetsBuf.contents(), targetOffsets, targetOffsets.count * 4)
        memcpy(targetLengthsBuf.contents(), targetLengths, targetLengths.count * 2)
        memcpy(matBuf.contents(), mat, mat.count)
        memcpy(paramsBuf.contents(), params, params.count * 2)
        memcpy(h0Buf.contents(), h0Values, h0Values.count * 2)
        memcpy(wsOffsetsBuf.contents(), wsOffsets, wsOffsets.count * 4)

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
        encoder.setBuffer(h0Buf, offset: 0, index: 8)
        encoder.setBuffer(resultsBuf, offset: 0, index: 9)
        encoder.setBuffer(wsBuf, offset: 0, index: 10)
        encoder.setBuffer(wsOffsetsBuf, offset: 0, index: 11)

        let threadgroupSize = min(64, pipeline.maxTotalThreadsPerThreadgroup)
        let gridSize = MTLSize(width: taskCount, height: 1, depth: 1)
        let tgSize = MTLSize(width: threadgroupSize, height: 1, depth: 1)
        encoder.dispatchThreads(gridSize, threadsPerThreadgroup: tgSize)
        encoder.endEncoding()

        cmdBuf.commit()
        cmdBuf.waitUntilCompleted()

        var results: [SWResult] = []
        results.reserveCapacity(taskCount)
        let resultPtr = resultsBuf.contents().bindMemory(to: Int32.self, capacity: taskCount * 6)

        for i in 0..<taskCount {
            let base = i * 6
            results.append(SWResult(
                score: resultPtr[base],
                queryEnd: resultPtr[base + 1],
                targetEnd: resultPtr[base + 2],
                globalTargetEnd: resultPtr[base + 3],
                globalScore: resultPtr[base + 4]
            ))
        }

        pool.release(allBuffers)
        return results
    }
}
#endif
