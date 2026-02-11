/// A Super-Maximal Exact Match (SMEM) interval in the FM-index.
///
/// Corresponds to `SMEM` struct in bwa-mem2's `FMI_search.h`.
/// Represents a BWT interval [k, l) with count s = l - k + 1 and query range [queryBegin, queryEnd).
public struct SMEM: Sendable, Equatable {
    /// Lower bound of the SA interval (inclusive)
    public var k: Int64
    /// Upper bound of the SA interval (inclusive)
    public var l: Int64
    /// Start position in the query (inclusive)
    public var queryBegin: Int32
    /// End position in the query (exclusive)
    public var queryEnd: Int32

    /// Number of occurrences (SA interval size)
    @inlinable
    public var count: Int64 { l - k + 1 }

    /// Length of the match in query coordinates
    @inlinable
    public var length: Int32 { queryEnd - queryBegin }

    public init(k: Int64, l: Int64, queryBegin: Int32, queryEnd: Int32) {
        self.k = k
        self.l = l
        self.queryBegin = queryBegin
        self.queryEnd = queryEnd
    }
}
