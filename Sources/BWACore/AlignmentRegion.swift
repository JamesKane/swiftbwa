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
