import Foundation
import BWACore

/// BWT data backed by an allocated checkpoint array.
public final class BWT: @unchecked Sendable {
    /// Checkpoint occurrence array
    public let checkpoints: UnsafeBufferPointer<CheckpointOCC>
    /// Cumulative base counts [A_before, C_before, G_before, T_before, total+1]
    /// After load, each count[i] has been incremented by 1 (matching bwa-mem2 load_index behavior)
    public let count: (Int64, Int64, Int64, Int64, Int64)
    /// Total length of the BWT (reference_seq_len in bwa-mem2, = 2*genome_len + 1)
    public let length: Int64
    /// Position of the sentinel character '$' in the BWT
    public let sentinelIndex: Int64

    /// One-hot mask array for rank queries (64 entries)
    /// one_hot_mask_array[0] = 0, [1] = 0x8000000000000000, [i] = [i-1] >> 1 | 0x8000000000000000
    public let oneHotMaskArray: [UInt64]

    // Owned allocation for cleanup
    private let ownedBase: UnsafeMutableRawPointer?

    public init(checkpoints: UnsafeBufferPointer<CheckpointOCC>,
                count: (Int64, Int64, Int64, Int64, Int64),
                length: Int64,
                sentinelIndex: Int64,
                ownedBase: UnsafeMutableRawPointer? = nil) {
        self.checkpoints = checkpoints
        self.count = count
        self.length = length
        self.sentinelIndex = sentinelIndex
        self.ownedBase = ownedBase

        // Build one-hot mask array (matches FMI_search.cpp:386-394)
        var masks = [UInt64](repeating: 0, count: 64)
        let high: UInt64 = 0x8000_0000_0000_0000
        masks[0] = 0
        if masks.count > 1 {
            masks[1] = high
            for i in 2..<64 {
                masks[i] = (masks[i - 1] >> 1) | high
            }
        }
        self.oneHotMaskArray = masks
    }

    deinit {
        if let base = ownedBase {
            base.deallocate()
        }
    }

    /// Get occurrence count of base `c` (0-3) at position `pos` in the BWT.
    /// Implements the GET_OCC macro from FMI_search.h:69-76.
    @inlinable
    public func occ(_ pos: Int64, _ c: Int) -> Int64 {
        let occID = Int(pos >> CP_SHIFT)
        let y = Int(pos & CP_MASK)

        let checkpoint = checkpoints[occID]
        let cpCount: Int64
        let bitstring: UInt64
        switch c {
        case 0: cpCount = checkpoint.counts.0; bitstring = checkpoint.bitstrings.0
        case 1: cpCount = checkpoint.counts.1; bitstring = checkpoint.bitstrings.1
        case 2: cpCount = checkpoint.counts.2; bitstring = checkpoint.bitstrings.2
        case 3: cpCount = checkpoint.counts.3; bitstring = checkpoint.bitstrings.3
        default: return 0
        }

        let matchMask = bitstring & oneHotMaskArray[y]
        return cpCount + Int64(matchMask.nonzeroBitCount)
    }
}
