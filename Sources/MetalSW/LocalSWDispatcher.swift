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
    /// Uses wavefront kernel (32-thread SIMD groups) when available, falls back to scalar.
    public static func dispatchBatch(
        tasks: [LocalSWTask],
        engine: MetalSWEngine
    ) async -> [LocalSWResult?] {
        guard !tasks.isEmpty else { return Array(repeating: nil, count: tasks.count) }

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
            fwdResults = await dispatchPassWavefront(tasks: tasks, pipeline: pipeline, engine: engine)
        } else {
            fwdResults = await dispatchPass(tasks: tasks, pipeline: pipeline, engine: engine)
        }

        // Build reverse tasks for start position recovery
        var revIndices: [Int] = []
        var revBatch: [LocalSWTask] = []
        for (i, task) in tasks.enumerated() {
            guard let fwd = fwdResults[i], fwd.score > 0 else { continue }
            let revQLen = Int(fwd.qEnd) + 1
            let revTLen = Int(fwd.tEnd) + 1
            var revQ = [UInt8](repeating: 0, count: revQLen)
            var revT = [UInt8](repeating: 0, count: revTLen)
            for k in 0..<revQLen { revQ[k] = task.query[revQLen - 1 - k] }
            for k in 0..<revTLen { revT[k] = task.target[revTLen - 1 - k] }
            revIndices.append(i)
            revBatch.append(LocalSWTask(query: revQ, target: revT, scoring: task.scoring))
        }

        // Reverse pass
        let revResults: [LocalSWPassResult?]
        if !revBatch.isEmpty {
            if useWavefront {
                revResults = await dispatchPassWavefront(tasks: revBatch, pipeline: pipeline, engine: engine)
            } else {
                revResults = await dispatchPass(tasks: revBatch, pipeline: pipeline, engine: engine)
            }
        } else {
            revResults = []
        }

        // Combine forward + reverse results
        var results: [LocalSWResult?] = Array(repeating: nil, count: tasks.count)
        for (revIdx, i) in revIndices.enumerated() {
            guard let fwd = fwdResults[i], let rev = revResults[revIdx] else { continue }
            results[i] = LocalSWResult(
                score: fwd.score,
                queryBegin: fwd.qEnd - rev.qEnd,
                queryEnd: fwd.qEnd,
                targetBegin: fwd.tEnd - rev.tEnd,
                targetEnd: fwd.tEnd
            )
        }
        return results
    }

    /// Read local SW pass results from GPU output buffer.
    private static func readPassResults(
        buffer: MTLBuffer, taskCount: Int
    ) -> [LocalSWPassResult?] {
        let ptr = buffer.contents().bindMemory(to: Int32.self, capacity: taskCount * 3)
        var results: [LocalSWPassResult?] = []
        results.reserveCapacity(taskCount)
        for i in 0..<taskCount {
            let base = i * 3
            let score = ptr[base]
            if score < 0 {
                results.append(nil)
            } else {
                results.append(LocalSWPassResult(score: score, qEnd: ptr[base + 1], tEnd: ptr[base + 2]))
            }
        }
        return results
    }

    /// Wavefront local SW pass: one SIMD group (32 threads) per task.
    private static func dispatchPassWavefront(
        tasks: [LocalSWTask],
        pipeline: MTLComputePipelineState,
        engine: MetalSWEngine
    ) async -> [LocalSWPassResult?] {
        let taskCount = tasks.count
        var packer = QueryTargetPacker()
        for task in tasks { packer.add(query: task.query, target: task.target) }

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let gapOE = scoring.gapOpenPenalty + scoring.gapExtendPenalty
        let params: [Int16] = [Int16(gapOE), Int16(scoring.gapExtendPenalty)]
        var numTasks: UInt32 = UInt32(taskCount)

        let pool = engine.bufferPool
        guard let qtBufs = packer.createBuffers(pool: pool),
              let extraBufs = acquireBuffers(pool: pool, sizes: [
                  mat.count, params.count * 2, taskCount * 3 * 4, 4
              ])
        else { return Array(repeating: nil, count: taskCount) }

        let allBuffers = qtBufs + extraBufs
        memcpy(extraBufs[0].contents(), mat, mat.count)
        memcpy(extraBufs[1].contents(), params, params.count * 2)
        memcpy(extraBufs[3].contents(), &numTasks, 4)

        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(allBuffers)
            return Array(repeating: nil, count: taskCount)
        }

        encoder.setComputePipelineState(pipeline)
        encoder.setBufferSequence(allBuffers)

        let simdGroupsPerTG = 4
        let threadsPerTG = simdGroupsPerTG * 32
        let numThreadgroups = (taskCount + simdGroupsPerTG - 1) / simdGroupsPerTG
        encoder.dispatchThreadgroups(
            MTLSize(width: numThreadgroups, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: threadsPerTG, height: 1, depth: 1)
        )
        encoder.endEncoding()

        await withCheckedContinuation { (continuation: CheckedContinuation<Void, Never>) in
            cmdBuf.addCompletedHandler { _ in continuation.resume() }
            cmdBuf.commit()
        }

        let results = readPassResults(buffer: extraBufs[2], taskCount: taskCount)
        pool.release(allBuffers)
        return results
    }

    /// Scalar single-direction local SW pass on GPU (fallback).
    private static func dispatchPass(
        tasks: [LocalSWTask],
        pipeline: MTLComputePipelineState,
        engine: MetalSWEngine
    ) async -> [LocalSWPassResult?] {
        let taskCount = tasks.count
        var packer = QueryTargetPacker()
        var wsOffsets: [UInt32] = []
        var totalWS: UInt32 = 0

        for task in tasks {
            packer.add(query: task.query, target: task.target)
            let qlen = task.query.count
            let wsBytes = UInt32(qlen + qlen + 5 * qlen)
            let aligned = (totalWS + 15) & ~15
            wsOffsets.append(aligned)
            totalWS = aligned + wsBytes
        }

        let scoring = tasks[0].scoring
        let mat = scoring.scoringMatrix()
        let gapOE = scoring.gapOpenPenalty + scoring.gapExtendPenalty
        let params: [Int16] = [Int16(gapOE), Int16(scoring.gapExtendPenalty), Int16(scoring.mismatchPenalty)]

        let pool = engine.bufferPool
        guard let qtBufs = packer.createBuffers(pool: pool),
              let extraBufs = acquireBuffers(pool: pool, sizes: [
                  mat.count, params.count * 2, taskCount * 3 * 4,
                  Int(totalWS), wsOffsets.count * 4
              ])
        else { return Array(repeating: nil, count: taskCount) }

        let allBuffers = qtBufs + extraBufs
        memcpy(extraBufs[0].contents(), mat, mat.count)
        memcpy(extraBufs[1].contents(), params, params.count * 2)
        memcpy(extraBufs[4].contents(), wsOffsets, wsOffsets.count * 4)

        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(allBuffers)
            return Array(repeating: nil, count: taskCount)
        }

        encoder.setComputePipelineState(pipeline)
        encoder.setBufferSequence(allBuffers)

        let tgSize = min(64, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(
            MTLSize(width: taskCount, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: tgSize, height: 1, depth: 1)
        )
        encoder.endEncoding()

        await withCheckedContinuation { (continuation: CheckedContinuation<Void, Never>) in
            cmdBuf.addCompletedHandler { _ in continuation.resume() }
            cmdBuf.commit()
        }

        let results = readPassResults(buffer: extraBufs[2], taskCount: taskCount)
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
