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
    /// Maps GPU dispatch position to original read index: sortedOrder[gpuIdx] = originalIdx.
    /// Nil when no reordering was applied.
    let sortedOrder: [Int]?
}

/// Batch dispatcher for GPU-accelerated SMEM forward extension.
public struct SMEMDispatcher: Sendable {

    /// Pack first 20 bases of a read into a UInt64 sort key (2 bits per base, N→0).
    /// Shorter reads are left-shifted so same-prefix reads sort together.
    @inline(__always)
    private static func sortKey(for query: [UInt8]) -> UInt64 {
        var key: UInt64 = 0
        let n = min(query.count, 20)
        for i in 0..<n {
            let base = query[i] < 4 ? UInt64(query[i]) : 0
            key = (key << 2) | base
        }
        key <<= UInt64((20 - n) * 2)
        return key
    }

    /// Dispatch GPU forward extension for a batch of reads.
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
        else { return nil }

        let taskCount = queries.count

        // Sort reads by prefix key for SIMD group cache coherence
        let sortedOrder: [Int]?
        let orderedQueries: [[UInt8]]
        if taskCount > 1 {
            let keys = queries.map { sortKey(for: $0) }
            let order = (0..<taskCount).sorted { keys[$0] < keys[$1] }
            sortedOrder = order
            orderedQueries = order.map { queries[$0] }
        } else {
            sortedOrder = nil
            orderedQueries = queries
        }

        var allQueries: [UInt8] = []
        var queryOffsets: [UInt32] = []
        var queryLengths: [UInt16] = []
        var resultOffsets: [UInt32] = []
        var totalPositions: Int = 0

        for query in orderedQueries {
            queryOffsets.append(UInt32(allQueries.count))
            queryLengths.append(UInt16(query.count))
            resultOffsets.append(UInt32(totalPositions))
            allQueries.append(contentsOf: query)
            totalPositions += query.count
        }

        let totalPos = max(totalPositions, 1)
        let pool = engine.bufferPool

        guard let bufs = acquireBuffers(pool: pool, sizes: [
            max(allQueries.count, 1), taskCount * 4, taskCount * 2,
            totalPos * 4, totalPos * 8, totalPos * 8,
            taskCount * 4, 4
        ]) else { return nil }

        if !allQueries.isEmpty { memcpy(bufs[0].contents(), allQueries, allQueries.count) }
        memcpy(bufs[1].contents(), queryOffsets, taskCount * 4)
        memcpy(bufs[2].contents(), queryLengths, taskCount * 2)
        memcpy(bufs[6].contents(), resultOffsets, taskCount * 4)
        var numTasks = UInt32(taskCount)
        memcpy(bufs[7].contents(), &numTasks, 4)

        guard let cmdBuf = engine.queue.makeCommandBuffer(),
              let encoder = cmdBuf.makeComputeCommandEncoder()
        else {
            pool.release(bufs)
            return nil
        }

        encoder.setComputePipelineState(pipeline)
        // Binding order: queries(0), offsets(1), lengths(2), bwt(3), counts(4), sentinel(5),
        // rightEnds(6→3 shift needed)... Actually the kernel bindings are specific:
        encoder.setBuffer(bufs[0], offset: 0, index: 0)   // queries
        encoder.setBuffer(bufs[1], offset: 0, index: 1)   // queryOffsets
        encoder.setBuffer(bufs[2], offset: 0, index: 2)   // queryLengths
        encoder.setBuffer(bwtBuf, offset: 0, index: 3)    // bwt checkpoints
        encoder.setBuffer(countsBuf, offset: 0, index: 4) // bwt counts
        encoder.setBuffer(sentBuf, offset: 0, index: 5)   // sentinel
        encoder.setBuffer(bufs[3], offset: 0, index: 6)   // rightEnds
        encoder.setBuffer(bufs[4], offset: 0, index: 7)   // kValues
        encoder.setBuffer(bufs[5], offset: 0, index: 8)   // sValues
        encoder.setBuffer(bufs[6], offset: 0, index: 9)   // resultOffsets
        encoder.setBuffer(bufs[7], offset: 0, index: 10)  // numTasks
        encoder.setBuffer(kmerBuf, offset: 0, index: 11)   // kmer hash

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
            rightEndsBuf: bufs[3], kValuesBuf: bufs[4], sValuesBuf: bufs[5],
            poolBuffers: bufs,
            queryLengths: queryLengths, resultOffsets: resultOffsets,
            totalPositions: totalPositions, taskCount: taskCount, pool: pool,
            sortedOrder: sortedOrder
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

        let results: [ForwardExtResult]

        if let sortedOrder = handle.sortedOrder {
            // Build inverse permutation and gather results in original order
            let count = handle.taskCount
            results = [ForwardExtResult](unsafeUninitializedCapacity: count) { buf, initCount in
                for gpuIdx in 0..<count {
                    let origIdx = sortedOrder[gpuIdx]
                    let offset = Int(handle.resultOffsets[gpuIdx])
                    let len = Int(handle.queryLengths[gpuIdx])
                    buf.initializeElement(at: origIdx, to: ForwardExtResult(
                        rightEnds: Array(UnsafeBufferPointer(start: rightEndPtr + offset, count: len)),
                        kValues: Array(UnsafeBufferPointer(start: kPtr + offset, count: len)),
                        sValues: Array(UnsafeBufferPointer(start: sPtr + offset, count: len))
                    ))
                }
                initCount = count
            }
        } else {
            var r: [ForwardExtResult] = []
            r.reserveCapacity(handle.taskCount)
            for i in 0..<handle.taskCount {
                let offset = Int(handle.resultOffsets[i])
                let len = Int(handle.queryLengths[i])
                r.append(ForwardExtResult(
                    rightEnds: Array(UnsafeBufferPointer(start: rightEndPtr + offset, count: len)),
                    kValues: Array(UnsafeBufferPointer(start: kPtr + offset, count: len)),
                    sValues: Array(UnsafeBufferPointer(start: sPtr + offset, count: len))
                ))
            }
            results = r
        }

        handle.pool.release(handle.poolBuffers)
        return results
    }

    /// Convert forward extension results to SMEMs.
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
                if rEnd <= Int32(i) { continue }
                if i > 0 && query[i - 1] < 4 && rEnd <= result.rightEnds[i - 1] { continue }
                let matchLen = rEnd - Int32(i)
                if matchLen < minSeedLen { continue }
                smems.append(SMEM(
                    k: result.kValues[i], l: result.kValues[i] + result.sValues[i] - 1,
                    queryBegin: Int32(i), queryEnd: rEnd
                ))
            }

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
