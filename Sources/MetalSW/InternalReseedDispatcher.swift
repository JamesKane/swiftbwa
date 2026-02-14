#if canImport(Metal)
import Metal
import BWACore

/// A single internal reseed task: find SMEMs at a midpoint with elevated minIntv.
public struct ReseedTask: Sendable {
    public let readIndex: Int
    public let startPos: UInt16
    public let minIntv: Int64
}

/// Handle for an in-flight GPU internal reseed dispatch.
public struct InternalReseedHandle: @unchecked Sendable {
    let cmdBuf: MTLCommandBuffer
    let outK: MTLBuffer
    let outL: MTLBuffer
    let outQBegin: MTLBuffer
    let outQEnd: MTLBuffer
    let outCount: MTLBuffer
    let poolBuffers: [MTLBuffer]
    let tasks: [ReseedTask]
    let pool: MetalBufferPool
}

/// Max SMEMs per task — must match MAX_OUT in internal_reseed.metal.
private let maxOutPerTask = 16

/// Batch dispatcher for GPU-accelerated internal reseeding.
/// Scans SMEMs for qualifying candidates, dispatches the internal_reseed
/// Metal kernel, and merges results back into per-read SMEM arrays.
public struct InternalReseedDispatcher: Sendable {

    /// Scan all reads' SMEMs and build reseed tasks for qualifying candidates.
    /// A candidate qualifies when: len >= splitLen AND count <= splitWidth.
    public static func buildTasks(
        allSMEMs: [[SMEM]],
        scoring: ScoringParameters
    ) -> [ReseedTask] {
        let splitLen = Int(Float(scoring.minSeedLength) * scoring.seedSplitRatio + 0.499)
        let splitWidth = Int64(scoring.splitWidth)
        var tasks: [ReseedTask] = []

        for (readIdx, smems) in allSMEMs.enumerated() {
            for smem in smems {
                let len = Int(smem.queryEnd - smem.queryBegin)
                if len < splitLen || smem.count > splitWidth { continue }
                let midpoint = UInt16((Int(smem.queryBegin) + Int(smem.queryEnd)) / 2)
                tasks.append(ReseedTask(
                    readIndex: readIdx,
                    startPos: midpoint,
                    minIntv: smem.count + 1
                ))
            }
        }
        return tasks
    }

    /// Dispatch GPU internal reseed kernel asynchronously.
    /// Returns a handle to collect results, or nil if dispatch cannot proceed.
    public static func dispatchBatchAsync(
        tasks: [ReseedTask],
        reads: [ReadSequence],
        scoring: ScoringParameters,
        engine: MetalSWEngine
    ) -> InternalReseedHandle? {
        guard let pipeline = engine.internalReseedPipeline,
              let bwtBuf = engine.bwtCheckpointBuffer,
              let countsBuf = engine.bwtCountsBuffer,
              let sentBuf = engine.bwtSentinelBuffer,
              !tasks.isEmpty
        else {
            return nil
        }

        let taskCount = tasks.count

        // Pack all reads' bases contiguously
        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []
        for read in reads {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(read.length))
            allQueries.append(contentsOf: read.bases)
        }

        // Pack per-task arrays
        var taskReadIdx: [UInt32] = []
        var taskStartPos: [UInt16] = []
        var taskMinIntv: [Int64] = []
        taskReadIdx.reserveCapacity(taskCount)
        taskStartPos.reserveCapacity(taskCount)
        taskMinIntv.reserveCapacity(taskCount)
        for task in tasks {
            taskReadIdx.append(UInt32(task.readIndex))
            taskStartPos.append(task.startPos)
            taskMinIntv.append(task.minIntv)
        }

        let pool = engine.bufferPool
        let readCount = reads.count

        // Acquire buffers
        guard let queriesBuf = pool.acquire(minSize: max(allQueries.count, 1)),
              let queryOffsetsBuf = pool.acquire(minSize: readCount * 4),
              let queryLengthsBuf = pool.acquire(minSize: readCount * 2),
              let taskReadIdxBuf = pool.acquire(minSize: taskCount * 4),
              let taskStartPosBuf = pool.acquire(minSize: taskCount * 2),
              let taskMinIntvBuf = pool.acquire(minSize: taskCount * 8),
              let outKBuf = pool.acquire(minSize: taskCount * maxOutPerTask * 8),
              let outLBuf = pool.acquire(minSize: taskCount * maxOutPerTask * 8),
              let outQBeginBuf = pool.acquire(minSize: taskCount * maxOutPerTask * 2),
              let outQEndBuf = pool.acquire(minSize: taskCount * maxOutPerTask * 2),
              let outCountBuf = pool.acquire(minSize: taskCount),
              let numTasksBuf = pool.acquire(minSize: 4),
              let minSeedLenBuf = pool.acquire(minSize: 2)
        else {
            return nil
        }

        let poolBuffers = [queriesBuf, queryOffsetsBuf, queryLengthsBuf,
                           taskReadIdxBuf, taskStartPosBuf, taskMinIntvBuf,
                           outKBuf, outLBuf, outQBeginBuf, outQEndBuf,
                           outCountBuf, numTasksBuf, minSeedLenBuf]

        // Copy input data
        if !allQueries.isEmpty {
            memcpy(queriesBuf.contents(), allQueries, allQueries.count)
        }
        memcpy(queryOffsetsBuf.contents(), queryOffsets, readCount * 4)
        memcpy(queryLengthsBuf.contents(), queryLengths, readCount * 2)
        memcpy(taskReadIdxBuf.contents(), taskReadIdx, taskCount * 4)
        memcpy(taskStartPosBuf.contents(), taskStartPos, taskCount * 2)
        memcpy(taskMinIntvBuf.contents(), taskMinIntv, taskCount * 8)
        var numTasks = UInt32(taskCount)
        memcpy(numTasksBuf.contents(), &numTasks, 4)
        var minSeedLen = Int16(scoring.minSeedLength)
        memcpy(minSeedLenBuf.contents(), &minSeedLen, 2)

        // Zero the output count buffer
        memset(outCountBuf.contents(), 0, taskCount)

        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(poolBuffers)
            return nil
        }

        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(queriesBuf, offset: 0, index: 0)
        encoder.setBuffer(queryOffsetsBuf, offset: 0, index: 1)
        encoder.setBuffer(queryLengthsBuf, offset: 0, index: 2)
        encoder.setBuffer(taskReadIdxBuf, offset: 0, index: 3)
        encoder.setBuffer(taskStartPosBuf, offset: 0, index: 4)
        encoder.setBuffer(taskMinIntvBuf, offset: 0, index: 5)
        encoder.setBuffer(bwtBuf, offset: 0, index: 6)
        encoder.setBuffer(countsBuf, offset: 0, index: 7)
        encoder.setBuffer(sentBuf, offset: 0, index: 8)
        encoder.setBuffer(outKBuf, offset: 0, index: 9)
        encoder.setBuffer(outLBuf, offset: 0, index: 10)
        encoder.setBuffer(outQBeginBuf, offset: 0, index: 11)
        encoder.setBuffer(outQEndBuf, offset: 0, index: 12)
        encoder.setBuffer(outCountBuf, offset: 0, index: 13)
        encoder.setBuffer(numTasksBuf, offset: 0, index: 14)
        encoder.setBuffer(minSeedLenBuf, offset: 0, index: 15)

        // One thread per task, 256 threads per threadgroup
        let threadsPerTG = min(256, pipeline.maxTotalThreadsPerThreadgroup)
        let numThreadgroups = (taskCount + threadsPerTG - 1) / threadsPerTG
        encoder.dispatchThreadgroups(
            MTLSize(width: numThreadgroups, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: threadsPerTG, height: 1, depth: 1)
        )
        encoder.endEncoding()
        cmdBuf.commit()

        return InternalReseedHandle(
            cmdBuf: cmdBuf,
            outK: outKBuf,
            outL: outLBuf,
            outQBegin: outQBeginBuf,
            outQEnd: outQEndBuf,
            outCount: outCountBuf,
            poolBuffers: poolBuffers,
            tasks: tasks,
            pool: pool
        )
    }

    /// Wait for GPU completion and merge results into per-read SMEM arrays.
    public static func collectAndMerge(
        handle: InternalReseedHandle,
        allSMEMs: inout [[SMEM]]
    ) {
        handle.cmdBuf.waitUntilCompleted()

        let taskCount = handle.tasks.count
        let countPtr = handle.outCount.contents().bindMemory(
            to: UInt8.self, capacity: taskCount)
        let kPtr = handle.outK.contents().bindMemory(
            to: Int64.self, capacity: taskCount * maxOutPerTask)
        let lPtr = handle.outL.contents().bindMemory(
            to: Int64.self, capacity: taskCount * maxOutPerTask)
        let qBeginPtr = handle.outQBegin.contents().bindMemory(
            to: Int16.self, capacity: taskCount * maxOutPerTask)
        let qEndPtr = handle.outQEnd.contents().bindMemory(
            to: Int16.self, capacity: taskCount * maxOutPerTask)

        // Track which reads got new SMEMs for re-sorting
        var affectedReads = Set<Int>()

        for tid in 0..<taskCount {
            let count = Int(countPtr[tid])
            if count == 0 { continue }

            let readIdx = handle.tasks[tid].readIndex
            let baseOff = tid * maxOutPerTask
            affectedReads.insert(readIdx)

            for i in 0..<count {
                let k = kPtr[baseOff + i]
                let l = lPtr[baseOff + i]
                let qBegin = Int32(qBeginPtr[baseOff + i])
                let qEnd = Int32(qEndPtr[baseOff + i])
                allSMEMs[readIdx].append(SMEM(
                    k: k, l: l, queryBegin: qBegin, queryEnd: qEnd
                ))
            }
        }

        // Re-sort affected reads' SMEM arrays
        for readIdx in affectedReads {
            allSMEMs[readIdx].sort {
                if $0.queryBegin != $1.queryBegin { return $0.queryBegin < $1.queryBegin }
                return $0.length > $1.length
            }
        }

        handle.pool.release(handle.poolBuffers)
    }
}
#endif
