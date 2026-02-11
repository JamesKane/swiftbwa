import BWACore

/// Checkpoint occurrence structure for FM-Index rank queries.
/// Each checkpoint covers 64 BWT characters (CP_BLOCK_SIZE=64).
/// Maps to `CP_OCC` in bwa-mem2's FMI_search.h.
public struct CheckpointOCC: Sendable {
    /// Cumulative counts of A, C, G, T up to this checkpoint
    public var counts: (Int64, Int64, Int64, Int64)
    /// One-hot encoded BWT bitstrings for 64 characters.
    /// Bit i is set if character at position (checkpoint_start + 63 - i) equals this base.
    public var bitstrings: (UInt64, UInt64, UInt64, UInt64)

    public init() {
        counts = (0, 0, 0, 0)
        bitstrings = (0, 0, 0, 0)
    }
}

// Constants
public let CP_BLOCK_SIZE: Int = 64
public let CP_SHIFT: Int = 6
public let CP_MASK: Int64 = 63
