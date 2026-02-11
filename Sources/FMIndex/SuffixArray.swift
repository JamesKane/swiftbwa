import Foundation
import BWACore

/// Compressed suffix array: stores SA entries as 8-bit high byte + 32-bit low word.
/// Supports both compressed (every 8th entry) and uncompressed modes.
public final class SuffixArray: @unchecked Sendable {
    public let msBytes: UnsafeBufferPointer<Int8>
    public let lsWords: UnsafeBufferPointer<UInt32>
    public let count: Int
    /// SA compression factor (0 = uncompressed, 3 = every 8th position)
    public let compressionShift: Int

    // Owned allocations for cleanup
    private let ownedMSBase: UnsafeMutableRawPointer?
    private let ownedLSBase: UnsafeMutableRawPointer?

    public init(msBytes: UnsafeBufferPointer<Int8>,
                lsWords: UnsafeBufferPointer<UInt32>,
                count: Int,
                compressionShift: Int = 0,
                ownedMSBase: UnsafeMutableRawPointer? = nil,
                ownedLSBase: UnsafeMutableRawPointer? = nil) {
        self.msBytes = msBytes
        self.lsWords = lsWords
        self.count = count
        self.compressionShift = compressionShift
        self.ownedMSBase = ownedMSBase
        self.ownedLSBase = ownedLSBase
    }

    deinit {
        if let base = ownedMSBase {
            base.deallocate()
        }
        if let base = ownedLSBase {
            base.deallocate()
        }
    }

    /// Get SA entry at a sampled position (direct array access).
    /// Only valid when `pos` is a multiple of `1 << compressionShift`.
    @inlinable
    public func entry(at pos: Int64) -> Int64 {
        let idx = Int(pos)
        let ms = Int64(msBytes[idx]) & 0xFF
        let ls = Int64(lsWords[idx])
        return (ms << 32) | ls
    }

    /// Resolve SA entry at an arbitrary BWT position using walk-backs.
    /// Implements bwa-mem2's `get_sa_entry64()`: walks through the BWT via
    /// LF mapping until hitting a sampled position, then adds the step count.
    @inlinable
    public func resolve(at pos: Int64, bwt: BWT) -> Int64 {
        let mask = Int64((1 << compressionShift) - 1)
        var p = pos
        var steps: Int64 = 0
        while (p & mask) != 0 {
            p = bwt.lfMapping(p)
            steps += 1
        }
        return entry(at: p >> Int64(compressionShift)) + steps
    }
}
