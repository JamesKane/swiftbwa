/// A chain of collinear seeds on the reference.
///
/// Corresponds to `mem_chain_t` in bwa-mem2's `bwamem.h`.
public struct MemChain: Sendable {
    /// Seeds in this chain, sorted by reference position
    public var seeds: [MemSeed]
    /// Chain weight (sum of seed lengths minus overlaps)
    public var weight: Int32
    /// Reference sequence ID
    public var rid: Int32
    /// Whether this chain is kept after filtering
    public var kept: Int32
    /// First seed index used for sorting
    public var first: Int32
    /// Reference start of the first seed
    public var pos: Int64
    /// Whether this chain maps to an ALT contig
    public var isAlt: Bool
    /// Fraction of read bases in highly repetitive seeds
    public var fracRep: Float = 0.0

    public init(seeds: [MemSeed] = [], weight: Int32 = 0, rid: Int32 = -1,
                kept: Int32 = 0, first: Int32 = -1, pos: Int64 = 0,
                isAlt: Bool = false) {
        self.seeds = seeds
        self.weight = weight
        self.rid = rid
        self.kept = kept
        self.first = first
        self.pos = pos
        self.isAlt = isAlt
    }
}
