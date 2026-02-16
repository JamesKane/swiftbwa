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

    /// Pack tasks into contiguous arrays and compute per-task workspace.
    private static func packTasks(
        _ tasks: [BandedSWTask], elementSize: Int
    ) -> (packer: QueryTargetPacker, h0Values: [Int16], wsOffsets: [UInt32], totalWS: UInt32) {
        var packer = QueryTargetPacker()
        var h0Values: [Int16] = []
        var wsOffsets: [UInt32] = []
        var totalWS: UInt32 = 0

        for task in tasks {
            packer.add(query: task.query, target: task.target)
            h0Values.append(Int16(clamping: task.h0))
            let qlen = task.query.count
            let wsBytes = UInt32((qlen + qlen + 5 * qlen) * elementSize)
            let aligned = (totalWS + 15) & ~15
            wsOffsets.append(aligned)
            totalWS = aligned + wsBytes
        }
        return (packer, h0Values, wsOffsets, totalWS)
    }

    /// Dispatch a batch of 8-bit banded SW tasks to the GPU.
    /// Returns SWResult for each task. Tasks that overflowed 8-bit have nil result.
    public static func dispatchBatch8(
        tasks: [BandedSWTask],
        engine: MetalSWEngine
    ) -> [SWResult?] {
        guard !tasks.isEmpty else { return [] }
        let taskCount = tasks.count
        let (packer, h0Values, wsOffsets, totalWS) = packTasks(tasks, elementSize: 1)

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let params: [Int16] = [
            Int16(scoring.gapOpenPenalty), Int16(scoring.gapExtendPenalty),
            Int16(scoring.gapOpenPenaltyDeletion), Int16(scoring.gapExtendPenaltyDeletion),
            Int16(scoring.zDrop), Int16(scoring.mismatchPenalty), Int16(tasks[0].w)
        ]

        let pool = engine.bufferPool
        guard let qtBufs = packer.createBuffers(pool: pool),
              let extraBufs = acquireBuffers(pool: pool, sizes: [
                  mat.count, params.count * 2, h0Values.count * 2,
                  taskCount * 6 * 4, Int(totalWS), wsOffsets.count * 4
              ])
        else { return Array(repeating: nil, count: taskCount) }

        let allBuffers = qtBufs + extraBufs
        memcpy(extraBufs[0].contents(), mat, mat.count)
        memcpy(extraBufs[1].contents(), params, params.count * 2)
        memcpy(extraBufs[2].contents(), h0Values, h0Values.count * 2)
        memcpy(extraBufs[5].contents(), wsOffsets, wsOffsets.count * 4)

        guard encodeAndDispatchSync(
            pipeline: engine.bandedSW8Pipeline, engine: engine,
            buffers: allBuffers, taskCount: taskCount
        ) else {
            pool.release(allBuffers)
            return Array(repeating: nil, count: taskCount)
        }

        let resultPtr = extraBufs[3].contents().bindMemory(to: Int32.self, capacity: taskCount * 6)
        var results: [SWResult?] = []
        results.reserveCapacity(taskCount)
        for i in 0..<taskCount {
            let base = i * 6
            let score = resultPtr[base]
            if score <= -2147483647 {
                results.append(nil)
            } else {
                results.append(SWResult(
                    score: score, queryEnd: resultPtr[base + 1],
                    targetEnd: resultPtr[base + 2], globalTargetEnd: resultPtr[base + 3],
                    globalScore: resultPtr[base + 4], maxOff: resultPtr[base + 5]
                ))
            }
        }
        pool.release(allBuffers)
        return results
    }

    /// Dispatch a batch of 16-bit banded SW tasks to the GPU.
    public static func dispatchBatch16(
        tasks: [BandedSWTask],
        engine: MetalSWEngine
    ) -> [SWResult] {
        guard !tasks.isEmpty, let pipeline = engine.bandedSW16Pipeline else { return [] }
        let taskCount = tasks.count
        let (packer, h0Values, wsOffsets, totalWS) = packTasks(tasks, elementSize: 2)

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let params: [Int16] = [
            Int16(scoring.gapOpenPenalty), Int16(scoring.gapExtendPenalty),
            Int16(scoring.gapOpenPenaltyDeletion), Int16(scoring.gapExtendPenaltyDeletion),
            Int16(scoring.zDrop), 0, Int16(tasks[0].w)
        ]

        let pool = engine.bufferPool
        guard let qtBufs = packer.createBuffers(pool: pool),
              let extraBufs = acquireBuffers(pool: pool, sizes: [
                  mat.count, params.count * 2, h0Values.count * 2,
                  taskCount * 6 * 4, Int(totalWS), wsOffsets.count * 4
              ])
        else { return [] }

        let allBuffers = qtBufs + extraBufs
        memcpy(extraBufs[0].contents(), mat, mat.count)
        memcpy(extraBufs[1].contents(), params, params.count * 2)
        memcpy(extraBufs[2].contents(), h0Values, h0Values.count * 2)
        memcpy(extraBufs[5].contents(), wsOffsets, wsOffsets.count * 4)

        guard encodeAndDispatchSync(
            pipeline: pipeline, engine: engine,
            buffers: allBuffers, taskCount: taskCount
        ) else {
            pool.release(allBuffers)
            return []
        }

        let resultPtr = extraBufs[3].contents().bindMemory(to: Int32.self, capacity: taskCount * 6)
        var results: [SWResult] = []
        results.reserveCapacity(taskCount)
        for i in 0..<taskCount {
            let base = i * 6
            results.append(SWResult(
                score: resultPtr[base], queryEnd: resultPtr[base + 1],
                targetEnd: resultPtr[base + 2], globalTargetEnd: resultPtr[base + 3],
                globalScore: resultPtr[base + 4], maxOff: resultPtr[base + 5]
            ))
        }
        pool.release(allBuffers)
        return results
    }
}
#endif
