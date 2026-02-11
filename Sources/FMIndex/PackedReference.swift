import Foundation
import BWACore

/// 2-bit packed reference sequence loaded from .pac file.
/// Each byte stores 4 bases: high bits first (bits 7-6 = first base, bits 1-0 = fourth base).
public final class PackedReference: @unchecked Sendable {
    /// Raw packed bytes
    public let data: UnsafeBufferPointer<UInt8>
    /// Total number of bases (not bytes)
    public let length: Int64

    // Owned allocation for cleanup
    private let ownedBase: UnsafeMutableRawPointer?

    public init(data: UnsafeBufferPointer<UInt8>, length: Int64,
                ownedBase: UnsafeMutableRawPointer? = nil) {
        self.data = data
        self.length = length
        self.ownedBase = ownedBase
    }

    deinit {
        if let base = ownedBase {
            base.deallocate()
        }
    }

    /// Get base at position (0-indexed). Returns 0-3 for A,C,G,T.
    @inlinable
    public func base(at pos: Int64) -> UInt8 {
        let byteIdx = Int(pos >> 2)
        let shift = (3 - Int(pos & 3)) << 1
        return (data[byteIdx] >> shift) & 3
    }

    /// Extract a subsequence as an array of 2-bit encoded bases.
    public func subsequence(from start: Int64, length: Int) -> [UInt8] {
        var result = [UInt8](repeating: 0, count: length)
        for i in 0..<length {
            result[i] = base(at: start + Int64(i))
        }
        return result
    }
}
