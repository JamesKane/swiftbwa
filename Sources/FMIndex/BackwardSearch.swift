import BWACore

/// FM-Index backward search operations.
/// Provides the core `backwardExt` operation that extends an SA interval by one character.
/// Reimplements bwa-mem2's `backwardExt()` from FMI_search.cpp:1025-1052.
public struct BackwardSearch: Sendable {

    /// Extend an SA interval backward by one character.
    ///
    /// Given a bidirectional interval (k, l, s) where:
    ///   - [k, k+s) is the SA interval in the forward BWT
    ///   - l tracks the reverse-direction SA interval
    /// compute the new interval for pattern cP (c prepended to P).
    ///
    /// This matches bwa-mem2's `backwardExt()` exactly:
    /// 1. Compute k[4] and s[4] for all 4 bases via Occ
    /// 2. Compute l[4] via a ladder from smem.l with sentinel adjustment
    /// 3. Return (k[a], l[a], s[a]) for the requested base
    @inlinable
    public static func backwardExt(
        bwt: BWT,
        interval: (k: Int64, l: Int64, s: Int64),
        base a: Int
    ) -> (k: Int64, l: Int64, s: Int64) {
        guard a < 4 else {
            return (0, 0, 0)
        }

        let sp = interval.k
        let ep = interval.k + interval.s

        // Compute k and s for all 4 bases
        var kAll: (Int64, Int64, Int64, Int64) = (0, 0, 0, 0)
        var sAll: (Int64, Int64, Int64, Int64) = (0, 0, 0, 0)

        // Base 0 (A)
        let occSP0 = bwt.occ(sp, 0)
        let occEP0 = bwt.occ(ep, 0)
        kAll.0 = bwt.count.0 + occSP0
        sAll.0 = occEP0 - occSP0

        // Base 1 (C)
        let occSP1 = bwt.occ(sp, 1)
        let occEP1 = bwt.occ(ep, 1)
        kAll.1 = bwt.count.1 + occSP1
        sAll.1 = occEP1 - occSP1

        // Base 2 (G)
        let occSP2 = bwt.occ(sp, 2)
        let occEP2 = bwt.occ(ep, 2)
        kAll.2 = bwt.count.2 + occSP2
        sAll.2 = occEP2 - occSP2

        // Base 3 (T)
        let occSP3 = bwt.occ(sp, 3)
        let occEP3 = bwt.occ(ep, 3)
        kAll.3 = bwt.count.3 + occSP3
        sAll.3 = occEP3 - occSP3

        // Compute l values via ladder from interval.l with sentinel adjustment
        let sentinelOffset: Int64 = (sp <= bwt.sentinelIndex && ep > bwt.sentinelIndex) ? 1 : 0
        var lAll: (Int64, Int64, Int64, Int64) = (0, 0, 0, 0)
        lAll.3 = interval.l + sentinelOffset
        lAll.2 = lAll.3 + sAll.3
        lAll.1 = lAll.2 + sAll.2
        lAll.0 = lAll.1 + sAll.1

        switch a {
        case 0: return (kAll.0, lAll.0, sAll.0)
        case 1: return (kAll.1, lAll.1, sAll.1)
        case 2: return (kAll.2, lAll.2, sAll.2)
        case 3: return (kAll.3, lAll.3, sAll.3)
        default: return (0, 0, 0)
        }
    }

    /// Initialize an SA interval for a single character.
    /// Returns the interval (k, l, s) for the bidirectional BWT.
    /// k = C[c], l = C[3-c], s = C[c+1] - C[c]
    @inlinable
    public static func initInterval(
        bwt: BWT,
        base c: Int
    ) -> (k: Int64, l: Int64, s: Int64) {
        guard c < 4 else { return (0, 0, 0) }

        let cCount = bwt.count(for: c)
        let cCountNext = bwt.count(forNext: c)
        let s = cCountNext - cCount
        return (cCount, cCount + s - 1, s)
    }
}
