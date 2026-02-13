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
    public let bufferPool: MetalBufferPool

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
            for name in ["banded_sw8", "banded_sw16", "local_sw", "local_sw_wavefront"] {
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
