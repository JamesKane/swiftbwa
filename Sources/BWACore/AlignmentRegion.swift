/// Hash function for deterministic tiebreaking among equal-score regions.
/// Matches bwa-mem2's `hash_64()` from utils.h.
public func hash64(_ key: UInt64) -> UInt64 {
    var k = key
    k &+= ~(k &<< 32)
    k ^= (k >> 22)
    k &+= ~(k &<< 13)
    k ^= (k >> 8)
    k &+= (k &<< 3)
    k ^= (k >> 15)
    k &+= ~(k &<< 27)
    k ^= (k >> 31)
    return k
}

/// An alignment region produced by Smith-Waterman extension of a chain.
///
/// Corresponds to `mem_alnreg_t` in bwa-mem2's `bwamem.h`.
public struct MemAlnReg: Sendable {
    /// Reference begin (in concatenated forward+reverse BWT coordinate space)
    public var rb: Int64
    /// Reference end
    public var re: Int64
    /// Query begin
    public var qb: Int32
    /// Query end
    public var qe: Int32
    /// Reference sequence ID
    public var rid: Int32
    /// Best SW score
    public var score: Int32
    /// Score from true scoring matrix (with N penalty)
    public var trueScore: Int32
    /// Best sub-optimal score
    public var sub: Int32
    /// SW score of a tandem hit (second-best local SW from mate rescue)
    public var csub: Int32 = 0
    /// Number of sub-optimal hits
    public var subN: Int32 = 0
    /// Alignment band width used
    public var w: Int32
    /// Number of seed bases covered
    public var seedCov: Int32
    /// Secondary alignment index (-1 if primary)
    public var secondary: Int32 = -1
    /// Hash for dedup
    public var hash: UInt64 = 0
    /// Number of seeds in the parent chain
    public var seedLen0: Int32 = 0
    /// Number of ambiguous bases
    public var nAmb: Int32 = 0
    /// Flap length (clipping)
    public var flapLen: Int32 = 0
    /// Whether this is an ALT hit
    public var isAlt: Bool = false
    /// Best ALT hit score (when this is a primary hit that has an ALT competitor)
    public var altSc: Int32 = 0
    /// Secondary index from Phase 1 (all-hits) marking, before primary-only re-ranking
    public var secondaryAll: Int32 = -1

    public init(rb: Int64 = 0, re: Int64 = 0, qb: Int32 = 0, qe: Int32 = 0,
                rid: Int32 = -1, score: Int32 = 0, trueScore: Int32 = 0,
                sub: Int32 = 0, w: Int32 = 0, seedCov: Int32 = 0) {
        self.rb = rb
        self.re = re
        self.qb = qb
        self.qe = qe
        self.rid = rid
        self.score = score
        self.trueScore = trueScore
        self.sub = sub
        self.w = w
        self.seedCov = seedCov
    }
}
