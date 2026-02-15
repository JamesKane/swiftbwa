#if canImport(Metal)
import Metal

/// Acquire multiple Metal buffers from a pool atomically.
/// If any acquisition fails, releases already-acquired buffers and returns nil.
func acquireBuffers(pool: MetalBufferPool, sizes: [Int]) -> [MTLBuffer]? {
    var buffers: [MTLBuffer] = []
    buffers.reserveCapacity(sizes.count)
    for size in sizes {
        guard let buf = pool.acquire(minSize: max(size, 1)) else {
            pool.release(buffers)
            return nil
        }
        buffers.append(buf)
    }
    return buffers
}

extension MTLComputeCommandEncoder {
    /// Set multiple buffers at sequential binding indices starting from 0.
    func setBufferSequence(_ buffers: [MTLBuffer]) {
        for (i, buf) in buffers.enumerated() {
            setBuffer(buf, offset: 0, index: i)
        }
    }
}

/// Packs variable-length query and target sequences into contiguous arrays
/// with offset/length metadata for GPU dispatch.
struct QueryTargetPacker {
    var allQueries: [UInt8] = []
    var queryOffsets: [UInt32] = []
    var queryLengths: [UInt16] = []
    var allTargets: [UInt8] = []
    var targetOffsets: [UInt32] = []
    var targetLengths: [UInt16] = []

    mutating func add(query: [UInt8], target: [UInt8]) {
        queryOffsets.append(UInt32(allQueries.count))
        queryLengths.append(UInt16(query.count))
        allQueries.append(contentsOf: query)
        targetOffsets.append(UInt32(allTargets.count))
        targetLengths.append(UInt16(target.count))
        allTargets.append(contentsOf: target)
    }

    /// Create and fill query+target Metal buffers (6 total).
    /// Copies data into the buffers. Returns nil if any acquisition fails.
    func createBuffers(pool: MetalBufferPool) -> [MTLBuffer]? {
        guard let bufs = acquireBuffers(pool: pool, sizes: [
            max(allQueries.count, 1), queryOffsets.count * 4, queryLengths.count * 2,
            max(allTargets.count, 1), targetOffsets.count * 4, targetLengths.count * 2
        ]) else { return nil }
        if !allQueries.isEmpty { memcpy(bufs[0].contents(), allQueries, allQueries.count) }
        memcpy(bufs[1].contents(), queryOffsets, queryOffsets.count * 4)
        memcpy(bufs[2].contents(), queryLengths, queryLengths.count * 2)
        if !allTargets.isEmpty { memcpy(bufs[3].contents(), allTargets, allTargets.count) }
        memcpy(bufs[4].contents(), targetOffsets, targetOffsets.count * 4)
        memcpy(bufs[5].contents(), targetLengths, targetLengths.count * 2)
        return bufs
    }
}

/// Create command buffer + encoder, set pipeline, bind all buffers, dispatch
/// with 1-thread-per-task grid, commit, and wait synchronously.
/// Returns false if command buffer creation fails (caller should release buffers).
func encodeAndDispatchSync(
    pipeline: MTLComputePipelineState,
    engine: MetalSWEngine,
    buffers: [MTLBuffer],
    taskCount: Int
) -> Bool {
    guard let cmdBuf = engine.queue.makeCommandBuffer(),
          let encoder = cmdBuf.makeComputeCommandEncoder()
    else { return false }

    encoder.setComputePipelineState(pipeline)
    encoder.setBufferSequence(buffers)

    let tgSize = min(64, pipeline.maxTotalThreadsPerThreadgroup)
    encoder.dispatchThreads(
        MTLSize(width: taskCount, height: 1, depth: 1),
        threadsPerThreadgroup: MTLSize(width: tgSize, height: 1, depth: 1)
    )
    encoder.endEncoding()
    cmdBuf.commit()
    cmdBuf.waitUntilCompleted()
    return true
}
#endif
