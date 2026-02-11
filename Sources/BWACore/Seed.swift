/// A positioned seed in the reference, created by looking up an SMEM in the suffix array.
///
/// Corresponds to `mem_seed_t` in bwa-mem2's `bwamem.h`.
public struct MemSeed: Sendable, Equatable {
    /// Start position in the reference (BWT coordinate)
    public var rbeg: Int64
    /// Start position in the query
    public var qbeg: Int32
    /// Length of the seed match
    public var len: Int32
    /// Seed score (typically = len * matchScore)
    public var score: Int32

    public init(rbeg: Int64, qbeg: Int32, len: Int32, score: Int32) {
        self.rbeg = rbeg
        self.qbeg = qbeg
        self.len = len
        self.score = score
    }
}
