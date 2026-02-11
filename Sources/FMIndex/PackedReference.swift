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

    /// Strand-aware sequence extraction in BWT coordinate space.
    ///
    /// Equivalent to bwa-mem2's `bns_fetch_seq()` / `bns_get_seq()`.
    ///
    /// - For forward strand (`bwtEnd <= genomeLength`): extracts directly.
    /// - For reverse strand (`bwtBegin >= genomeLength`): converts to forward coordinates,
    ///   extracts, reverses, and complements each base.
    /// - Clamps to chromosome boundaries using `metadata.annotations[rid]`.
    ///
    /// - Parameters:
    ///   - bwtBegin: Start position in BWT coordinate space `[0, 2*genomeLength)`
    ///   - bwtEnd: End position (exclusive) in BWT coordinate space
    ///   - genomeLength: Forward genome length (half of total BWT length)
    ///   - metadata: Reference metadata for chromosome boundary lookup
    /// - Returns: Extracted bases and clamped coordinates, or nil if region is empty after clamping
    public func fetchSequence(
        bwtBegin: Int64, bwtEnd: Int64,
        genomeLength: Int64, metadata: ReferenceMetadata
    ) -> (bases: [UInt8], clampedBegin: Int64, clampedEnd: Int64, rid: Int32)? {
        guard bwtBegin < bwtEnd else { return nil }

        let isReverse = bwtBegin >= genomeLength

        if !isReverse {
            // Forward strand: clamp to chromosome boundary
            let rid = metadata.sequenceID(for: bwtBegin)
            let ann = metadata.annotations[Int(rid)]
            let chromEnd = ann.offset + Int64(ann.length)

            let cb = max(bwtBegin, ann.offset)
            let ce = min(bwtEnd, chromEnd)
            guard cb < ce else { return nil }

            let len = Int(ce - cb)
            let bases = subsequence(from: cb, length: len)
            return (bases, cb, ce, rid)
        } else {
            // Reverse strand: convert to forward coordinates
            // BWT reverse pos p maps to forward pos (2*genomeLength - 1 - p)
            let fwdBegin = 2 * genomeLength - bwtEnd
            let fwdEnd = 2 * genomeLength - bwtBegin

            let rid = metadata.sequenceID(for: fwdBegin)
            let ann = metadata.annotations[Int(rid)]
            let chromEnd = ann.offset + Int64(ann.length)

            let cb = max(fwdBegin, ann.offset)
            let ce = min(fwdEnd, chromEnd)
            guard cb < ce else { return nil }

            let len = Int(ce - cb)
            var bases = subsequence(from: cb, length: len)

            // Reverse and complement
            bases.reverse()
            for i in 0..<bases.count {
                bases[i] = 3 - bases[i]  // A↔T, C↔G
            }

            // Convert clamped forward coords back to BWT reverse coords
            let clampedBwtBegin = 2 * genomeLength - ce
            let clampedBwtEnd = 2 * genomeLength - cb
            return (bases, clampedBwtBegin, clampedBwtEnd, rid)
        }
    }
}
