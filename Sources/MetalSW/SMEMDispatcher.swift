#if canImport(Metal)
import Metal
import BWACore

/// Result of GPU forward extension for a single read.
/// Contains per-position (rightEnd, k, s) arrays.
public struct ForwardExtResult: Sendable {
    public let rightEnds: [Int32]    // length = query length
    public let kValues: [Int64]
    public let sValues: [Int64]

    public init(rightEnds: [Int32], kValues: [Int64], sValues: [Int64]) {
        self.rightEnds = rightEnds
        self.kValues = kValues
        self.sValues = sValues
    }
}

/// Handle for an in-flight GPU SMEM dispatch. Holds the command buffer and
/// result buffers so the caller can overlap CPU work before collecting results.
public struct SMEMDispatchHandle: @unchecked Sendable {
    let cmdBuf: MTLCommandBuffer
    let rightEndsBuf: MTLBuffer
    let kValuesBuf: MTLBuffer
    let sValuesBuf: MTLBuffer
    let poolBuffers: [MTLBuffer]
    let queryLengths: [UInt16]
    let resultOffsets: [UInt32]
    let totalPositions: Int
    let taskCount: Int
    let pool: MetalBufferPool
}

/// Batch dispatcher for GPU-accelerated SMEM forward extension.
/// Dispatches the smem_forward Metal kernel and converts results to SMEMs.
public struct SMEMDispatcher: Sendable {

    /// Dispatch GPU forward extension for a batch of reads.
    /// Returns one ForwardExtResult per read.
    public static func dispatchBatch(
        queries: [[UInt8]],
        engine: MetalSWEngine
    ) -> [ForwardExtResult] {
        guard let handle = dispatchBatchAsync(queries: queries, engine: engine) else {
            return queries.map { q in
                ForwardExtResult(
                    rightEnds: [Int32](repeating: 0, count: q.count),
                    kValues: [Int64](repeating: 0, count: q.count),
                    sValues: [Int64](repeating: 0, count: q.count)
                )
            }
        }
        return collectResults(handle: handle)
    }

    /// Dispatch GPU forward extension without waiting for completion.
    /// Returns a handle that must be passed to `collectResults` to retrieve data.
    /// Returns nil if the dispatch could not be set up (missing pipeline, empty queries, etc).
    public static func dispatchBatchAsync(
        queries: [[UInt8]],
        engine: MetalSWEngine
    ) -> SMEMDispatchHandle? {
        guard let pipeline = engine.smemForwardPipeline,
              let bwtBuf = engine.bwtCheckpointBuffer,
              let countsBuf = engine.bwtCountsBuffer,
              let sentBuf = engine.bwtSentinelBuffer,
              let kmerBuf = engine.kmerHashBuffer,
              !queries.isEmpty
        else {
            return nil
        }

        let taskCount = queries.count

        // Pack queries contiguously
        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []
        var resultOffsets: [UInt32] = []
        var totalPositions: Int = 0

        for query in queries {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(query.count))
            resultOffsets.append(UInt32(totalPositions))
            allQueries.append(contentsOf: query)
            totalPositions += query.count
        }

        // Pad to at least 1 to avoid zero-length buffer issues
        let totalPos = max(totalPositions, 1)

        let pool = engine.bufferPool

        guard let queriesBuf = pool.acquire(minSize: max(allQueries.count, 1)),
              let queryOffsetsBuf = pool.acquire(minSize: taskCount * 4),
              let queryLengthsBuf = pool.acquire(minSize: taskCount * 2),
              let rightEndsBuf = pool.acquire(minSize: totalPos * 4),
              let kValuesBuf = pool.acquire(minSize: totalPos * 8),
              let sValuesBuf = pool.acquire(minSize: totalPos * 8),
              let resultOffsetsBuf = pool.acquire(minSize: taskCount * 4),
              let numTasksBuf = pool.acquire(minSize: 4)
        else {
            return nil
        }

        let poolBuffers = [queriesBuf, queryOffsetsBuf, queryLengthsBuf,
                           rightEndsBuf, kValuesBuf, sValuesBuf,
                           resultOffsetsBuf, numTasksBuf]

        // Copy input data
        if !allQueries.isEmpty {
            memcpy(queriesBuf.contents(), allQueries, allQueries.count)
        }
        memcpy(queryOffsetsBuf.contents(), queryOffsets, taskCount * 4)
        memcpy(queryLengthsBuf.contents(), queryLengths, taskCount * 2)
        memcpy(resultOffsetsBuf.contents(), resultOffsets, taskCount * 4)
        var numTasks = UInt32(taskCount)
        memcpy(numTasksBuf.contents(), &numTasks, 4)

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
        encoder.setBuffer(bwtBuf, offset: 0, index: 3)
        encoder.setBuffer(countsBuf, offset: 0, index: 4)
        encoder.setBuffer(sentBuf, offset: 0, index: 5)
        encoder.setBuffer(rightEndsBuf, offset: 0, index: 6)
        encoder.setBuffer(kValuesBuf, offset: 0, index: 7)
        encoder.setBuffer(sValuesBuf, offset: 0, index: 8)
        encoder.setBuffer(resultOffsetsBuf, offset: 0, index: 9)
        encoder.setBuffer(numTasksBuf, offset: 0, index: 10)
        encoder.setBuffer(kmerBuf, offset: 0, index: 11)

        // One SIMD group (32 threads) per read, 4 SIMD groups per threadgroup
        let simdGroupsPerTG = 4
        let threadsPerTG = simdGroupsPerTG * 32
        let numThreadgroups = (taskCount + simdGroupsPerTG - 1) / simdGroupsPerTG

        encoder.dispatchThreadgroups(
            MTLSize(width: numThreadgroups, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: threadsPerTG, height: 1, depth: 1)
        )
        encoder.endEncoding()

        cmdBuf.commit()

        return SMEMDispatchHandle(
            cmdBuf: cmdBuf,
            rightEndsBuf: rightEndsBuf,
            kValuesBuf: kValuesBuf,
            sValuesBuf: sValuesBuf,
            poolBuffers: poolBuffers,
            queryLengths: queryLengths,
            resultOffsets: resultOffsets,
            totalPositions: totalPositions,
            taskCount: taskCount,
            pool: pool
        )
    }

    /// Wait for an in-flight GPU dispatch to complete and read back results.
    public static func collectResults(handle: SMEMDispatchHandle) -> [ForwardExtResult] {
        handle.cmdBuf.waitUntilCompleted()

        let rightEndPtr = handle.rightEndsBuf.contents().bindMemory(
            to: Int32.self, capacity: handle.totalPositions)
        let kPtr = handle.kValuesBuf.contents().bindMemory(
            to: Int64.self, capacity: handle.totalPositions)
        let sPtr = handle.sValuesBuf.contents().bindMemory(
            to: Int64.self, capacity: handle.totalPositions)

        var results: [ForwardExtResult] = []
        results.reserveCapacity(handle.taskCount)

        for i in 0..<handle.taskCount {
            let offset = Int(handle.resultOffsets[i])
            let len = Int(handle.queryLengths[i])

            let rightEnds = Array(UnsafeBufferPointer(start: rightEndPtr + offset, count: len))
            let kVals = Array(UnsafeBufferPointer(start: kPtr + offset, count: len))
            let sVals = Array(UnsafeBufferPointer(start: sPtr + offset, count: len))

            results.append(ForwardExtResult(
                rightEnds: rightEnds,
                kValues: kVals,
                sValues: sVals
            ))
        }

        handle.pool.release(handle.poolBuffers)
        return results
    }

    /// Convert forward extension results to SMEMs.
    /// Identifies left-maximal positions where rightEnd[i] > rightEnd[i-1],
    /// and applies minSeedLen filter.
    public static func extractSMEMs(
        from results: [ForwardExtResult],
        queries: [[UInt8]],
        minSeedLen: Int32
    ) -> [[SMEM]] {
        var allSMEMs: [[SMEM]] = []
        allSMEMs.reserveCapacity(results.count)

        for (readIdx, result) in results.enumerated() {
            let query = queries[readIdx]
            let readLen = query.count
            var smems: [SMEM] = []

            for i in 0..<readLen {
                let rEnd = result.rightEnds[i]

                // No match at this position (N base or empty)
                if rEnd <= Int32(i) { continue }

                // Left-maximality check: match is left-maximal if:
                // - i == 0 (leftmost position), or
                // - previous base is N (can't extend left), or
                // - rightEnd[i] > rightEnd[i-1] (this match extends further right)
                if i > 0 && query[i - 1] < 4 && rEnd <= result.rightEnds[i - 1] {
                    continue
                }

                let matchLen = rEnd - Int32(i)
                if matchLen < minSeedLen { continue }

                let k = result.kValues[i]
                let s = result.sValues[i]

                smems.append(SMEM(
                    k: k,
                    l: k + s - 1,
                    queryBegin: Int32(i),
                    queryEnd: rEnd
                ))
            }

            // Sort by query begin, then length descending (matching CPU convention)
            smems.sort {
                if $0.queryBegin != $1.queryBegin { return $0.queryBegin < $1.queryBegin }
                return $0.length > $1.length
            }

            allSMEMs.append(smems)
        }

        return allSMEMs
    }
}
#endif
