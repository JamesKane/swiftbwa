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
public struct InternalReseedDispatcher: Sendable {

    /// Scan all reads' SMEMs and build reseed tasks for qualifying candidates.
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
                tasks.append(ReseedTask(readIndex: readIdx, startPos: midpoint, minIntv: smem.count + 1))
            }
        }
        return tasks
    }

    /// Dispatch GPU internal reseed kernel asynchronously.
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
        else { return nil }

        let taskCount = tasks.count
        let readCount = reads.count

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
        guard let bufs = acquireBuffers(pool: pool, sizes: [
            max(allQueries.count, 1), readCount * 4, readCount * 2,
            taskCount * 4, taskCount * 2, taskCount * 8,
            taskCount * maxOutPerTask * 8, taskCount * maxOutPerTask * 8,
            taskCount * maxOutPerTask * 2, taskCount * maxOutPerTask * 2,
            taskCount, 4, 2
        ]) else { return nil }

        // Copy input data
        if !allQueries.isEmpty { memcpy(bufs[0].contents(), allQueries, allQueries.count) }
        memcpy(bufs[1].contents(), queryOffsets, readCount * 4)
        memcpy(bufs[2].contents(), queryLengths, readCount * 2)
        memcpy(bufs[3].contents(), taskReadIdx, taskCount * 4)
        memcpy(bufs[4].contents(), taskStartPos, taskCount * 2)
        memcpy(bufs[5].contents(), taskMinIntv, taskCount * 8)
        var numTasks = UInt32(taskCount)
        memcpy(bufs[11].contents(), &numTasks, 4)
        var minSeedLen = Int16(scoring.minSeedLength)
        memcpy(bufs[12].contents(), &minSeedLen, 2)
        memset(bufs[10].contents(), 0, taskCount)  // zero output count

        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(bufs)
            return nil
        }

        encoder.setComputePipelineState(pipeline)
        // Bind: queries(0), qOffsets(1), qLengths(2), taskReadIdx(3), taskStartPos(4), taskMinIntv(5)
        for i in 0..<6 { encoder.setBuffer(bufs[i], offset: 0, index: i) }
        // BWT data at indices 6-8
        encoder.setBuffer(bwtBuf, offset: 0, index: 6)
        encoder.setBuffer(countsBuf, offset: 0, index: 7)
        encoder.setBuffer(sentBuf, offset: 0, index: 8)
        // Output buffers at indices 9-13, numTasks at 14, minSeedLen at 15
        for i in 6..<13 { encoder.setBuffer(bufs[i], offset: 0, index: i + 3) }

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
            outK: bufs[6], outL: bufs[7], outQBegin: bufs[8], outQEnd: bufs[9], outCount: bufs[10],
            poolBuffers: bufs, tasks: tasks, pool: pool
        )
    }

    /// Wait for GPU completion and merge results into per-read SMEM arrays.
    public static func collectAndMerge(
        handle: InternalReseedHandle,
        allSMEMs: inout [[SMEM]]
    ) {
        handle.cmdBuf.waitUntilCompleted()

        let taskCount = handle.tasks.count
        let countPtr = handle.outCount.contents().bindMemory(to: UInt8.self, capacity: taskCount)
        let kPtr = handle.outK.contents().bindMemory(to: Int64.self, capacity: taskCount * maxOutPerTask)
        let lPtr = handle.outL.contents().bindMemory(to: Int64.self, capacity: taskCount * maxOutPerTask)
        let qBeginPtr = handle.outQBegin.contents().bindMemory(to: Int16.self, capacity: taskCount * maxOutPerTask)
        let qEndPtr = handle.outQEnd.contents().bindMemory(to: Int16.self, capacity: taskCount * maxOutPerTask)

        var affectedReads = Set<Int>()
        for tid in 0..<taskCount {
            let count = Int(countPtr[tid])
            if count == 0 { continue }
            let readIdx = handle.tasks[tid].readIndex
            let baseOff = tid * maxOutPerTask
            affectedReads.insert(readIdx)
            for i in 0..<count {
                allSMEMs[readIdx].append(SMEM(
                    k: kPtr[baseOff + i], l: lPtr[baseOff + i],
                    queryBegin: Int32(qBeginPtr[baseOff + i]), queryEnd: Int32(qEndPtr[baseOff + i])
                ))
            }
        }

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
