import BWACore

/// FM-Index backward search operations.
/// Provides the core `backwardExt` operation that extends an SA interval by one character.
/// Reimplements bwa-mem2's `backwardExt()` from FMI_search.cpp.
public struct BackwardSearch: Sendable {

    /// Extend an SA interval backward by one character.
    ///
    /// Given an interval [k, k+s) in the suffix array representing all suffixes prefixed
    /// by pattern P, compute the new interval for pattern cP (c prepended to P).
    ///
    /// This is the core FM-index operation: uses the C[] array and Occ() function.
    /// New interval: k' = C[c] + Occ(c, k-1), l' = C[c] + Occ(c, k+s-1) - 1
    ///
    /// Matches bwa-mem2's `backwardExt()` private method in FMI_search.
    @inlinable
    public static func backwardExt(
        bwt: BWT,
        interval: (k: Int64, l: Int64, s: Int64),
        base c: Int
    ) -> (k: Int64, l: Int64, s: Int64) {
        guard c < 4 else {
            // N base: interval collapses
            return (0, 0, 0)
        }

        let k = interval.k
        let s = interval.s

        // C[c] = count[c] (cumulative count of bases < c, plus 1 from load adjustment)
        let cCount: Int64
        switch c {
        case 0: cCount = bwt.count.0
        case 1: cCount = bwt.count.1
        case 2: cCount = bwt.count.2
        case 3: cCount = bwt.count.3
        default: cCount = 0
        }

        // Handle sentinel: if position falls on or crosses sentinel, adjust
        var occK: Int64
        var occKS: Int64

        if k <= bwt.sentinelIndex {
            occK = k > 0 ? bwt.occ(k - 1, c) : 0
        } else {
            occK = k > 0 ? bwt.occ(k - 1, c) : 0
        }

        let ks = k + s - 1
        if ks <= bwt.sentinelIndex {
            occKS = bwt.occ(ks, c)
        } else {
            occKS = bwt.occ(ks, c)
        }

        let newK = cCount + occK
        let newS = occKS - occK
        let newL = newK + newS - 1

        return (newK, newL, newS)
    }

    /// Initialize an SA interval for a single character.
    /// Returns the interval [C[c], C[c+1]) representing all suffixes starting with c.
    @inlinable
    public static func initInterval(
        bwt: BWT,
        base c: Int
    ) -> (k: Int64, l: Int64, s: Int64) {
        guard c < 4 else { return (0, 0, 0) }

        let cCount: Int64
        let cCountNext: Int64
        switch c {
        case 0: cCount = bwt.count.0; cCountNext = bwt.count.1
        case 1: cCount = bwt.count.1; cCountNext = bwt.count.2
        case 2: cCount = bwt.count.2; cCountNext = bwt.count.3
        case 3: cCount = bwt.count.3; cCountNext = bwt.count.4
        default: cCount = 0; cCountNext = 0
        }

        let s = cCountNext - cCount
        return (cCount, cCount + s - 1, s)
    }
}
