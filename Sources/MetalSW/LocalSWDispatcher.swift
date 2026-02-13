#if canImport(Metal)
import Metal
import BWACore

/// A local SW task to dispatch to the GPU.
public struct LocalSWTask: Sendable {
    public let query: [UInt8]
    public let target: [UInt8]
    public let scoring: ScoringParameters

    public init(query: [UInt8], target: [UInt8], scoring: ScoringParameters) {
        self.query = query
        self.target = target
        self.scoring = scoring
    }
}

/// Result of a local SW GPU dispatch (single pass).
public struct LocalSWPassResult: Sendable {
    public var score: Int32
    public var qEnd: Int32
    public var tEnd: Int32
}

/// Batch dispatcher for local Smith-Waterman on Metal GPU.
/// Handles the two-pass (forward + reverse) approach for start position recovery.
public struct LocalSWDispatcher: Sendable {

    /// Perform full local SW with start position recovery for a batch of tasks.
    /// Returns LocalSWResult for each task, nil if overflow or no alignment found.
    /// Uses wavefront kernel (32-thread SIMD groups) when available, falls back to scalar.
    public static func dispatchBatch(
        tasks: [LocalSWTask],
        engine: MetalSWEngine
    ) -> [LocalSWResult?] {
        guard !tasks.isEmpty else {
            return Array(repeating: nil, count: tasks.count)
        }

        // Prefer wavefront kernel (32-way parallelism per alignment)
        let useWavefront = engine.localSWWavefrontPipeline != nil
        let pipeline: MTLComputePipelineState
        if useWavefront, let wf = engine.localSWWavefrontPipeline {
            pipeline = wf
        } else if let scalar = engine.localSWPipeline {
            pipeline = scalar
        } else {
            return Array(repeating: nil, count: tasks.count)
        }

        // Forward pass
        let fwdResults: [LocalSWPassResult?]
        if useWavefront {
            fwdResults = dispatchPassWavefront(tasks: tasks, pipeline: pipeline, engine: engine)
        } else {
            fwdResults = dispatchPass(tasks: tasks, pipeline: pipeline, engine: engine)
        }

        // Build reverse tasks for start position recovery
        var revTasks: [LocalSWTask?] = []
        for (i, task) in tasks.enumerated() {
            guard let fwd = fwdResults[i], fwd.score > 0 else {
                revTasks.append(nil)
                continue
            }
            // Reverse prefix of query up to qEnd+1 and target up to tEnd+1
            let revQLen = Int(fwd.qEnd) + 1
            let revTLen = Int(fwd.tEnd) + 1
            var revQ = [UInt8](repeating: 0, count: revQLen)
            var revT = [UInt8](repeating: 0, count: revTLen)
            for k in 0..<revQLen { revQ[k] = task.query[revQLen - 1 - k] }
            for k in 0..<revTLen { revT[k] = task.target[revTLen - 1 - k] }
            revTasks.append(LocalSWTask(query: revQ, target: revT, scoring: task.scoring))
        }

        // Collect non-nil reverse tasks
        var revIndices: [Int] = []
        var revBatch: [LocalSWTask] = []
        for (i, task) in revTasks.enumerated() {
            if let t = task {
                revIndices.append(i)
                revBatch.append(t)
            }
        }

        // Reverse pass
        let revResults: [LocalSWPassResult?]
        if !revBatch.isEmpty {
            if useWavefront {
                revResults = dispatchPassWavefront(tasks: revBatch, pipeline: pipeline, engine: engine)
            } else {
                revResults = dispatchPass(tasks: revBatch, pipeline: pipeline, engine: engine)
            }
        } else {
            revResults = []
        }

        // Combine forward + reverse results
        var results: [LocalSWResult?] = Array(repeating: nil, count: tasks.count)
        var revIdx = 0
        for i in revIndices {
            guard let fwd = fwdResults[i], let rev = revResults[revIdx] else {
                revIdx += 1
                continue
            }
            results[i] = LocalSWResult(
                score: fwd.score,
                queryBegin: fwd.qEnd - rev.qEnd,
                queryEnd: fwd.qEnd,
                targetBegin: fwd.tEnd - rev.tEnd,
                targetEnd: fwd.tEnd
            )
            revIdx += 1
        }

        return results
    }

    /// Wavefront local SW pass: one SIMD group (32 threads) per task.
    private static func dispatchPassWavefront(
        tasks: [LocalSWTask],
        pipeline: MTLComputePipelineState,
        engine: MetalSWEngine
    ) -> [LocalSWPassResult?] {
        let taskCount = tasks.count

        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []
        var allTargets: [UInt8] = []
        var targetOffsets: [UInt32] = []
        var targetLengths: [UInt16] = []

        for task in tasks {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(task.query.count))
            allQueries.append(contentsOf: task.query)

            targetOffsets.append(UInt32(allTargets.count))
            targetLengths.append(UInt16(task.target.count))
            allTargets.append(contentsOf: task.target)
        }

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let gapOE = scoring.gapOpenPenalty + scoring.gapExtendPenalty
        let gapE = scoring.gapExtendPenalty
        let params: [Int16] = [Int16(gapOE), Int16(gapE)]
        var numTasks: UInt32 = UInt32(taskCount)

        let pool = engine.bufferPool

        guard let queriesBuf = pool.acquire(minSize: max(allQueries.count, 1)),
              let queryOffsetsBuf = pool.acquire(minSize: queryOffsets.count * 4),
              let queryLengthsBuf = pool.acquire(minSize: queryLengths.count * 2),
              let targetsBuf = pool.acquire(minSize: max(allTargets.count, 1)),
              let targetOffsetsBuf = pool.acquire(minSize: targetOffsets.count * 4),
              let targetLengthsBuf = pool.acquire(minSize: targetLengths.count * 2),
              let matBuf = pool.acquire(minSize: mat.count),
              let paramsBuf = pool.acquire(minSize: params.count * 2),
              let resultsBuf = pool.acquire(minSize: taskCount * 3 * 4),
              let numTasksBuf = pool.acquire(minSize: 4)
        else {
            return Array(repeating: nil, count: taskCount)
        }

        let allBuffers = [queriesBuf, queryOffsetsBuf, queryLengthsBuf,
                          targetsBuf, targetOffsetsBuf, targetLengthsBuf,
                          matBuf, paramsBuf, resultsBuf, numTasksBuf]

        if !allQueries.isEmpty { memcpy(queriesBuf.contents(), allQueries, allQueries.count) }
        memcpy(queryOffsetsBuf.contents(), queryOffsets, queryOffsets.count * 4)
        memcpy(queryLengthsBuf.contents(), queryLengths, queryLengths.count * 2)
        if !allTargets.isEmpty { memcpy(targetsBuf.contents(), allTargets, allTargets.count) }
        memcpy(targetOffsetsBuf.contents(), targetOffsets, targetOffsets.count * 4)
        memcpy(targetLengthsBuf.contents(), targetLengths, targetLengths.count * 2)
        memcpy(matBuf.contents(), mat, mat.count)
        memcpy(paramsBuf.contents(), params, params.count * 2)
        memcpy(numTasksBuf.contents(), &numTasks, 4)

        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(allBuffers)
            return Array(repeating: nil, count: taskCount)
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
        encoder.setBuffer(resultsBuf, offset: 0, index: 8)
        encoder.setBuffer(numTasksBuf, offset: 0, index: 9)

        // One SIMD group (32 threads) per task, 4 SIMD groups per threadgroup
        let simdGroupsPerTG = 4
        let threadsPerTG = simdGroupsPerTG * 32
        let numThreadgroups = (taskCount + simdGroupsPerTG - 1) / simdGroupsPerTG

        encoder.dispatchThreadgroups(
            MTLSize(width: numThreadgroups, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: threadsPerTG, height: 1, depth: 1)
        )
        encoder.endEncoding()

        cmdBuf.commit()
        cmdBuf.waitUntilCompleted()

        var results: [LocalSWPassResult?] = []
        results.reserveCapacity(taskCount)
        let resultPtr = resultsBuf.contents().bindMemory(to: Int32.self, capacity: taskCount * 3)

        for i in 0..<taskCount {
            let base = i * 3
            let score = resultPtr[base]
            if score < 0 {
                results.append(nil)
            } else {
                results.append(LocalSWPassResult(
                    score: score,
                    qEnd: resultPtr[base + 1],
                    tEnd: resultPtr[base + 2]
                ))
            }
        }

        pool.release(allBuffers)
        return results
    }

    /// Scalar single-direction local SW pass on GPU (fallback).
    private static func dispatchPass(
        tasks: [LocalSWTask],
        pipeline: MTLComputePipelineState,
        engine: MetalSWEngine
    ) -> [LocalSWPassResult?] {
        let taskCount = tasks.count

        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []
        var allTargets: [UInt8] = []
        var targetOffsets: [UInt32] = []
        var targetLengths: [UInt16] = []
        var wsOffsets: [UInt32] = []
        var totalWorkspaceBytes: UInt32 = 0

        for task in tasks {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(task.query.count))
            allQueries.append(contentsOf: task.query)

            targetOffsets.append(UInt32(allTargets.count))
            targetLengths.append(UInt16(task.target.count))
            allTargets.append(contentsOf: task.target)

            let qlen = task.query.count
            let wsBytes = UInt32(qlen + qlen + 5 * qlen)  // H + E + profile

            let aligned = (totalWorkspaceBytes + 15) & ~15
            wsOffsets.append(aligned)
            totalWorkspaceBytes = aligned + wsBytes
        }

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let gapOE = scoring.gapOpenPenalty + scoring.gapExtendPenalty
        let gapE = scoring.gapExtendPenalty
        let bias = scoring.mismatchPenalty
        let params: [Int16] = [Int16(gapOE), Int16(gapE), Int16(bias)]

        let pool = engine.bufferPool

        guard let queriesBuf = pool.acquire(minSize: max(allQueries.count, 1)),
              let queryOffsetsBuf = pool.acquire(minSize: queryOffsets.count * 4),
              let queryLengthsBuf = pool.acquire(minSize: queryLengths.count * 2),
              let targetsBuf = pool.acquire(minSize: max(allTargets.count, 1)),
              let targetOffsetsBuf = pool.acquire(minSize: targetOffsets.count * 4),
              let targetLengthsBuf = pool.acquire(minSize: targetLengths.count * 2),
              let matBuf = pool.acquire(minSize: mat.count),
              let paramsBuf = pool.acquire(minSize: params.count * 2),
              let resultsBuf = pool.acquire(minSize: taskCount * 3 * 4),
              let wsBuf = pool.acquire(minSize: max(Int(totalWorkspaceBytes), 1)),
              let wsOffsetsBuf = pool.acquire(minSize: wsOffsets.count * 4)
        else {
            return Array(repeating: nil, count: taskCount)
        }

        let allBuffers = [queriesBuf, queryOffsetsBuf, queryLengthsBuf,
                          targetsBuf, targetOffsetsBuf, targetLengthsBuf,
                          matBuf, paramsBuf, resultsBuf, wsBuf, wsOffsetsBuf]

        if !allQueries.isEmpty { memcpy(queriesBuf.contents(), allQueries, allQueries.count) }
        memcpy(queryOffsetsBuf.contents(), queryOffsets, queryOffsets.count * 4)
        memcpy(queryLengthsBuf.contents(), queryLengths, queryLengths.count * 2)
        if !allTargets.isEmpty { memcpy(targetsBuf.contents(), allTargets, allTargets.count) }
        memcpy(targetOffsetsBuf.contents(), targetOffsets, targetOffsets.count * 4)
        memcpy(targetLengthsBuf.contents(), targetLengths, targetLengths.count * 2)
        memcpy(matBuf.contents(), mat, mat.count)
        memcpy(paramsBuf.contents(), params, params.count * 2)
        memcpy(wsOffsetsBuf.contents(), wsOffsets, wsOffsets.count * 4)

        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(allBuffers)
            return Array(repeating: nil, count: taskCount)
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
        encoder.setBuffer(resultsBuf, offset: 0, index: 8)
        encoder.setBuffer(wsBuf, offset: 0, index: 9)
        encoder.setBuffer(wsOffsetsBuf, offset: 0, index: 10)

        let threadgroupSize = min(64, pipeline.maxTotalThreadsPerThreadgroup)
        let gridSize = MTLSize(width: taskCount, height: 1, depth: 1)
        let tgSize = MTLSize(width: threadgroupSize, height: 1, depth: 1)
        encoder.dispatchThreads(gridSize, threadsPerThreadgroup: tgSize)
        encoder.endEncoding()

        cmdBuf.commit()
        cmdBuf.waitUntilCompleted()

        var results: [LocalSWPassResult?] = []
        results.reserveCapacity(taskCount)
        let resultPtr = resultsBuf.contents().bindMemory(to: Int32.self, capacity: taskCount * 3)

        for i in 0..<taskCount {
            let base = i * 3
            let score = resultPtr[base]
            if score < 0 {
                results.append(nil)  // overflow
            } else {
                results.append(LocalSWPassResult(
                    score: score,
                    qEnd: resultPtr[base + 1],
                    tEnd: resultPtr[base + 2]
                ))
            }
        }

        pool.release(allBuffers)
        return results
    }
}

/// Re-export LocalSWResult so MetalSW doesn't need to depend on Alignment.
public struct LocalSWResult: Sendable {
    public var score: Int32
    public var queryBegin: Int32
    public var queryEnd: Int32
    public var targetBegin: Int32
    public var targetEnd: Int32

    public init(score: Int32, queryBegin: Int32, queryEnd: Int32,
                targetBegin: Int32, targetEnd: Int32) {
        self.score = score
        self.queryBegin = queryBegin
        self.queryEnd = queryEnd
        self.targetBegin = targetBegin
        self.targetEnd = targetEnd
    }
}
#endif
