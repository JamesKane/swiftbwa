import BWACore

/// Result of global banded alignment (ksw_global2 equivalent).
public struct GlobalAlignResult: Sendable {
    /// Alignment score
    public var score: Int32
    /// Packed CIGAR operations (length << 4 | op)
    public var cigar: [UInt32]
    /// Edit distance (number of mismatches + indel bases)
    public var nm: Int32

    public init(score: Int32 = 0, cigar: [UInt32] = [], nm: Int32 = 0) {
        self.score = score
        self.cigar = cigar
        self.nm = nm
    }
}

/// Result of CIGAR generation with position adjustment and soft-clipping.
public struct CIGARResult: Sendable {
    /// Final packed CIGAR (including soft-clips)
    public var cigar: [UInt32]
    /// Edit distance
    public var nm: Int32
    /// Alignment score from global DP
    public var score: Int32
    /// Adjusted reference position (after leading deletion squeeze)
    public var pos: Int64

    public init(cigar: [UInt32] = [], nm: Int32 = 0, score: Int32 = 0, pos: Int64 = 0) {
        self.cigar = cigar
        self.nm = nm
        self.score = score
        self.pos = pos
    }
}
