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
    public static func dispatchBatch(
        tasks: [GlobalSWTask],
        engine: MetalSWEngine
    ) -> [GlobalSWResult] {
        guard !tasks.isEmpty, let pipeline = engine.globalSWPipeline else { return [] }
        let taskCount = tasks.count

        var packer = QueryTargetPacker()
        var trueScoreValues: [Int16] = []
        var initialWValues: [Int16] = []
        var wsOffsets: [UInt32] = []
        var totalWS: UInt32 = 0

        for task in tasks {
            packer.add(query: task.query, target: task.target)
            trueScoreValues.append(Int16(clamping: task.trueScore))
            initialWValues.append(Int16(clamping: task.initialW))

            let qlen = task.query.count
            let tlen = task.target.count
            let maxNCol = qlen
            let wsBytes = UInt32(5 * qlen * 4 + (qlen + 1) * 4 + (qlen + 1) * 4 + tlen * maxNCol)
            let aligned = (totalWS + 15) & ~15
            wsOffsets.append(aligned)
            totalWS = aligned + wsBytes
        }

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let params: [Int16] = [
            Int16(scoring.gapOpenPenalty), Int16(scoring.gapExtendPenalty),
            Int16(scoring.gapOpenPenaltyDeletion), Int16(scoring.gapExtendPenaltyDeletion)
        ]

        let pool = engine.bufferPool
        guard let qtBufs = packer.createBuffers(pool: pool),
              let extraBufs = acquireBuffers(pool: pool, sizes: [
                  mat.count, params.count * 2, trueScoreValues.count * 2,
                  initialWValues.count * 2, taskCount * MAX_CIGAR_OPS * 4,
                  taskCount * 4, taskCount * 4, Int(totalWS), wsOffsets.count * 4
              ])
        else { return [] }

        let allBuffers = qtBufs + extraBufs
        memcpy(extraBufs[0].contents(), mat, mat.count)
        memcpy(extraBufs[1].contents(), params, params.count * 2)
        memcpy(extraBufs[2].contents(), trueScoreValues, trueScoreValues.count * 2)
        memcpy(extraBufs[3].contents(), initialWValues, initialWValues.count * 2)
        memcpy(extraBufs[8].contents(), wsOffsets, wsOffsets.count * 4)

        guard encodeAndDispatchSync(
            pipeline: pipeline, engine: engine,
            buffers: allBuffers, taskCount: taskCount
        ) else {
            pool.release(allBuffers)
            return []
        }

        // Read results
        let cigarPtr = extraBufs[4].contents().bindMemory(to: Int32.self, capacity: taskCount * MAX_CIGAR_OPS)
        let cigarLenPtr = extraBufs[5].contents().bindMemory(to: Int32.self, capacity: taskCount)
        let scorePtr = extraBufs[6].contents().bindMemory(to: Int32.self, capacity: taskCount)

        var results: [GlobalSWResult] = []
        results.reserveCapacity(taskCount)
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
