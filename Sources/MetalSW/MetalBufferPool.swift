#if canImport(Metal)
import Metal
import Foundation

/// Reusable MTLBuffer pool to avoid per-dispatch allocation overhead.
/// Thread-safe via os_unfair_lock.
public final class MetalBufferPool: @unchecked Sendable {
    private let device: MTLDevice
    private var pool: [MTLBuffer] = []
    private let lock = NSLock()

    /// Maximum number of buffers to keep in the pool.
    private let maxPoolSize = 32

    init(device: MTLDevice) {
        self.device = device
    }

    /// Acquire a buffer with at least `minSize` bytes.
    /// Reuses an existing buffer if one is large enough, otherwise allocates new.
    public func acquire(minSize: Int) -> MTLBuffer? {
        lock.lock()
        defer { lock.unlock() }

        // Find smallest buffer that fits
        var bestIdx: Int? = nil
        var bestSize = Int.max
        for i in 0..<pool.count {
            let bufSize = pool[i].length
            if bufSize >= minSize && bufSize < bestSize {
                bestIdx = i
                bestSize = bufSize
            }
        }

        if let idx = bestIdx {
            let buffer = pool.remove(at: idx)
            return buffer
        }

        // Allocate new — use shared storage for zero-copy on Apple Silicon
        return device.makeBuffer(length: max(minSize, 4096), options: .storageModeShared)
    }

    /// Return a buffer to the pool for reuse.
    public func release(_ buffer: MTLBuffer) {
        lock.lock()
        defer { lock.unlock() }

        if pool.count < maxPoolSize {
            pool.append(buffer)
        }
        // Otherwise let it deallocate
    }

    /// Release multiple buffers back to the pool.
    public func release(_ buffers: [MTLBuffer]) {
        lock.lock()
        defer { lock.unlock() }

        for buffer in buffers {
            if pool.count < maxPoolSize {
                pool.append(buffer)
            }
        }
    }
}
#endif
