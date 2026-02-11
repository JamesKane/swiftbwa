/// Scoring parameters matching BWA-MEM defaults from `mem_opt_init()`.
public struct ScoringParameters: Sendable {
    public var matchScore: Int32 = 1
    public var mismatchPenalty: Int32 = 4
    public var gapOpenPenalty: Int32 = 6
    public var gapExtendPenalty: Int32 = 1
    /// Separate gap open/extend for deletions (default: same as above)
    public var gapOpenPenaltyDeletion: Int32 = 6
    public var gapExtendPenaltyDeletion: Int32 = 1
    public var bandWidth: Int32 = 100
    public var zDrop: Int32 = 100
    public var minSeedLength: Int32 = 19
    public var maxOccurrences: Int32 = 500
    public var minOutputScore: Int32 = 30
    public var seedSplitRatio: Float = 1.5
    public var maxChainGap: Int32 = 10000
    public var maskLevel: Float = 0.50
    public var unpairedPenalty: Int32 = 17
    public var maxMatesw: Int32 = 50
    public var splitWidth: Int32 = 10
    public var maxInsert: Int32 = 10000
    /// Chunk size for batch processing
    public var chunkSize: Int = 10_000_000
    /// Number of threads
    public var numThreads: Int = 1
    /// Penalty clipping
    public var penClip5: Int32 = 5
    public var penClip3: Int32 = 5
    /// Flag bits (controls algorithm behavior)
    public var flag: Int32 = 0

    /// -M: Mark shorter split hits as secondary instead of supplementary
    public static let flagNoMulti: Int32 = 0x10
    /// -Y: Use soft clipping for supplementary alignments
    public static let flagSoftClip: Int32 = 0x200

    public init() {}

    /// Build the scoring matrix (5x5: A,C,G,T,N)
    public func scoringMatrix() -> [Int8] {
        var mat = [Int8](repeating: 0, count: 25)
        for i in 0..<4 {
            for j in 0..<4 {
                mat[i * 5 + j] = i == j ? Int8(matchScore) : Int8(-mismatchPenalty)
            }
            mat[i * 5 + 4] = -1  // vs N
            mat[4 * 5 + i] = -1  // N vs
        }
        mat[24] = -1  // N vs N
        return mat
    }
}
