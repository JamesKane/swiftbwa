#if canImport(Metal)
import Metal
import Foundation

/// Singleton Metal compute engine for Smith-Waterman GPU acceleration.
/// Returns nil on systems without Metal support (Linux, CI).
public final class MetalSWEngine: @unchecked Sendable {
    public static let shared: MetalSWEngine? = {
        do {
            return try MetalSWEngine()
        } catch {
            #if DEBUG
            fputs("[MetalSW] Engine init failed: \(error)\n", stderr)
            #endif
            return nil
        }
    }()

    public let device: MTLDevice
    public let queue: MTLCommandQueue
    let bandedSW8Pipeline: MTLComputePipelineState
    var bandedSW16Pipeline: MTLComputePipelineState?
    var localSWPipeline: MTLComputePipelineState?
    var localSWWavefrontPipeline: MTLComputePipelineState?
    public private(set) var smemForwardPipeline: MTLComputePipelineState?
    public private(set) var internalReseedPipeline: MTLComputePipelineState?
    public private(set) var globalSWPipeline: MTLComputePipelineState?
    public let bufferPool: MetalBufferPool

    // Cached BWT buffers for GPU seeding (created once via setupBWT, reused)
    public private(set) var bwtCheckpointBuffer: MTLBuffer?
    public private(set) var bwtCountsBuffer: MTLBuffer?
    public private(set) var bwtSentinelBuffer: MTLBuffer?

    /// Pre-computed 12-mer hash table for GPU seeding (384MB, built once).
    public private(set) var kmerHashBuffer: MTLBuffer?

    /// Minimum batch size to justify GPU dispatch overhead.
    public static let minBatchSize = 32

    private init() throws {
        guard let device = MTLCreateSystemDefaultDevice() else {
            throw MetalSWError.noDevice
        }
        guard let queue = device.makeCommandQueue() else {
            throw MetalSWError.noCommandQueue
        }
        self.device = device
        self.queue = queue
        self.bufferPool = MetalBufferPool(device: device)

        // Load Metal library: try pre-compiled metallib first, then compile from source
        let library: MTLLibrary
        if let bundleURL = Bundle.module.url(forResource: "default", withExtension: "metallib"),
           let lib = try? device.makeLibrary(URL: bundleURL) {
            library = lib
        } else {
            // Compile from .metal source files bundled via SPM .process("Kernels")
            // SPM flattens the resource directory, so files are at bundle root
            var sources: [String] = []
            for name in ["banded_sw8", "banded_sw16", "local_sw", "local_sw_wavefront", "smem_forward", "internal_reseed", "global_sw"] {
                if let url = Bundle.module.url(forResource: name, withExtension: "metal") {
                    sources.append(try String(contentsOf: url, encoding: .utf8))
                }
            }
            guard !sources.isEmpty else {
                throw MetalSWError.libraryNotFound
            }
            let combined = sources.joined(separator: "\n")
            library = try device.makeLibrary(source: combined, options: nil)
        }

        guard let sw8Fn = library.makeFunction(name: "banded_sw8") else {
            throw MetalSWError.functionNotFound("banded_sw8")
        }
        self.bandedSW8Pipeline = try device.makeComputePipelineState(function: sw8Fn)

        // Optional pipelines — loaded when available
        if let sw16Fn = library.makeFunction(name: "banded_sw16") {
            self.bandedSW16Pipeline = try device.makeComputePipelineState(function: sw16Fn)
        }
        if let localFn = library.makeFunction(name: "local_sw") {
            self.localSWPipeline = try device.makeComputePipelineState(function: localFn)
        }
        if let wfFn = library.makeFunction(name: "local_sw_wavefront") {
            self.localSWWavefrontPipeline = try device.makeComputePipelineState(function: wfFn)
        }
        if let smemFn = library.makeFunction(name: "smem_forward") {
            self.smemForwardPipeline = try device.makeComputePipelineState(function: smemFn)
        }
        if let reseedFn = library.makeFunction(name: "internal_reseed") {
            self.internalReseedPipeline = try device.makeComputePipelineState(function: reseedFn)
        }
        if let globalFn = library.makeFunction(name: "global_sw") {
            self.globalSWPipeline = try device.makeComputePipelineState(function: globalFn)
        }
    }

    /// Set up BWT checkpoint buffer for GPU seeding. Called once at pipeline start.
    /// Creates Metal buffers wrapping/copying the BWT checkpoint array, counts, and sentinel.
    ///
    /// - Parameters:
    ///   - checkpointPointer: Raw pointer to the CheckpointOCC array
    ///   - checkpointByteCount: Total size in bytes of the checkpoint array
    ///   - counts: The 5 cumulative base counts (A, C, G, T, total+1)
    ///   - sentinelIndex: Position of the '$' sentinel in the BWT
    public func setupBWT(
        checkpointPointer: UnsafeRawPointer,
        checkpointByteCount: Int,
        counts: (Int64, Int64, Int64, Int64, Int64),
        sentinelIndex: Int64
    ) {
        guard bwtCheckpointBuffer == nil else { return }  // already set up

        // Try zero-copy buffer if pointer is page-aligned
        let pageSize = Int(getpagesize())
        let pointerValue = Int(bitPattern: checkpointPointer)
        if pointerValue % pageSize == 0 {
            // Page-aligned: try bytesNoCopy for zero-copy GPU access
            let mutablePtr = UnsafeMutableRawPointer(mutating: checkpointPointer)
            bwtCheckpointBuffer = device.makeBuffer(
                bytesNoCopy: mutablePtr,
                length: checkpointByteCount,
                options: .storageModeShared,
                deallocator: nil  // BWT owns the memory
            )
        }

        // Fallback: copy the checkpoint data
        if bwtCheckpointBuffer == nil {
            bwtCheckpointBuffer = device.makeBuffer(
                bytes: checkpointPointer,
                length: checkpointByteCount,
                options: .storageModeShared
            )
        }

        // Counts buffer: 5 × Int64 = 40 bytes
        var countsArray: [Int64] = [counts.0, counts.1, counts.2, counts.3, counts.4]
        bwtCountsBuffer = device.makeBuffer(
            bytes: &countsArray,
            length: 5 * MemoryLayout<Int64>.size,
            options: .storageModeShared
        )

        // Sentinel buffer: 1 × Int64 = 8 bytes
        var sentinel = sentinelIndex
        bwtSentinelBuffer = device.makeBuffer(
            bytes: &sentinel,
            length: MemoryLayout<Int64>.size,
            options: .storageModeShared
        )
    }

    /// Set up the 12-mer hash table for GPU seeding.
    /// Takes a raw pointer to 4^12 = 16,777,216 entries of (k, l, s) × Int64 = 24 bytes each.
    /// Total: 384MB, built once by the caller and cached here.
    public func setupKmerHash(pointer: UnsafeRawPointer, byteCount: Int) {
        guard kmerHashBuffer == nil else { return }
        kmerHashBuffer = device.makeBuffer(
            bytes: pointer,
            length: byteCount,
            options: .storageModeShared
        )
    }
}

public enum MetalSWError: Error {
    case noDevice
    case noCommandQueue
    case libraryNotFound
    case functionNotFound(String)
    case encodingFailed
    case commandBufferFailed(String)
}
#endif
