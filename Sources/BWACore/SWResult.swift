/// Result of a banded Smith-Waterman extension.
public struct SWResult: Sendable, Equatable {
    /// Best alignment score
    public var score: Int32
    /// Query end position (length of query consumed)
    public var queryEnd: Int32
    /// Target (reference) end position (length of target consumed)
    public var targetEnd: Int32
    /// Global target end (when reaching end of query)
    public var globalTargetEnd: Int32
    /// Global score (alignment score when reaching end of query)
    public var globalScore: Int32

    public init(score: Int32 = 0, queryEnd: Int32 = 0, targetEnd: Int32 = 0,
                globalTargetEnd: Int32 = 0, globalScore: Int32 = 0) {
        self.score = score
        self.queryEnd = queryEnd
        self.targetEnd = targetEnd
        self.globalTargetEnd = globalTargetEnd
        self.globalScore = globalScore
    }
}
