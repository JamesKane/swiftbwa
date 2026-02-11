import BWACore

/// Finds Super-Maximal Exact Matches (SMEMs) using backward search on the FM-Index.
///
/// Reimplements `getSMEMsOnePosOneThread()` and `getSMEMsAllPosOneThread()` from
/// bwa-mem2's FMI_search.cpp (lines 496-724).
///
/// An SMEM is a maximal exact match that cannot be extended in either direction without
/// losing all occurrences. The algorithm works by:
/// 1. Forward extension (rightward) from each query position
/// 2. Backward extension (leftward) to find all maximal matches
public struct SMEMFinder: Sendable {

    /// Find all SMEMs for a query sequence at all positions.
    ///
    /// Reimplements `getSMEMsAllPosOneThread()` (FMI_search.cpp:672-724).
    /// Iterates over query positions, calling `findSMEMsAtPosition` until all done.
    ///
    /// - Parameters:
    ///   - query: 2-bit encoded query bases (A=0, C=1, G=2, T=3, N=4)
    ///   - bwt: The BWT structure for backward search
    ///   - minSeedLen: Minimum seed length to report
    /// - Returns: Array of SMEMs sorted by query position
    public static func findAllSMEMs(
        query: [UInt8],
        bwt: BWT,
        minSeedLen: Int32
    ) -> [SMEM] {
        let readLength = query.count
        guard readLength > 0 else { return [] }

        var allSMEMs: [SMEM] = []
        var pos: Int = 0

        while pos < readLength {
            let (smems, nextPos) = findSMEMsAtPosition(
                query: query,
                bwt: bwt,
                startPos: pos,
                minSeedLen: minSeedLen,
                minIntv: 1
            )
            allSMEMs.append(contentsOf: smems)
            pos = nextPos
        }

        // Sort by query begin position, then by length descending
        allSMEMs.sort {
            if $0.queryBegin != $1.queryBegin {
                return $0.queryBegin < $1.queryBegin
            }
            return $0.length > $1.length
        }

        return allSMEMs
    }

    /// Find SMEMs starting from a specific query position.
    ///
    /// Reimplements `getSMEMsOnePosOneThread()` (FMI_search.cpp:496-670).
    ///
    /// The algorithm:
    /// 1. Forward phase: extend rightward from `startPos`, tracking SA interval changes.
    ///    Each time the SA interval shrinks, record the previous interval as a candidate.
    /// 2. Backward phase: for each candidate, extend leftward, emitting SMEMs when the
    ///    interval drops below threshold.
    ///
    /// - Parameters:
    ///   - query: 2-bit encoded query
    ///   - bwt: BWT structure
    ///   - startPos: Position to start from
    ///   - minSeedLen: Minimum seed length
    ///   - minIntv: Minimum SA interval size (occurrence threshold)
    /// - Returns: (smems, nextPosition) tuple
    public static func findSMEMsAtPosition(
        query: [UInt8],
        bwt: BWT,
        startPos: Int,
        minSeedLen: Int32,
        minIntv: Int64
    ) -> (smems: [SMEM], nextPos: Int) {
        let readLength = query.count
        var smems: [SMEM] = []
        var nextPos = startPos + 1

        let a = query[startPos]
        guard a < 4 else {
            // N base at start position - skip
            return (smems, nextPos)
        }

        // Initialize SA interval for the first character
        // Forward search uses reverse complement: count[a] for k, count[3-a] for l
        var curK = bwt.count(for: Int(a))
        var curL = bwt.count(for: Int(3 - a))
        var curS = bwt.count(forNext: Int(a)) - bwt.count(for: Int(a))

        // Track previous intervals (candidates for backward extension)
        var prevArray: [(k: Int64, l: Int64, s: Int64, m: Int32, n: Int32)] = []
        var numPrev = 0

        // === Forward phase: extend rightward ===
        var j = startPos + 1
        while j < readLength {
            let base = query[j]
            nextPos = j + 1

            guard base < 4 else { break }

            // Forward extension is backward extension on reverse complement BWT
            // Swap k <-> l, use complement base (3-a)
            let swapped = (k: curL, l: curK, s: curS)
            let ext = BackwardSearch.backwardExt(bwt: bwt, interval: swapped, base: Int(3 - base))
            // Swap back
            let newK = ext.l
            let newL = ext.k
            let newS = ext.s

            // If interval size changed, record previous as candidate
            if newS != curS {
                prevArray.append((curK, curL, curS, Int32(startPos), Int32(j - 1)))
                numPrev += 1
            }

            // If new interval is too small, stop forward extension
            if newS < minIntv {
                nextPos = j
                break
            }

            curK = newK
            curL = newL
            curS = newS
            j += 1
        }

        // Record the final forward interval if it meets threshold
        if curS >= minIntv {
            prevArray.append((curK, curL, curS, Int32(startPos), Int32(j - 1)))
            numPrev += 1
        }

        // Reverse prevArray so longest match is first (for backward extension priority)
        prevArray.reverse()

        // === Backward phase: extend leftward from each candidate ===
        var prev = prevArray
        var currentJ = readLength  // nolint: used for tracking backward extension position
        _ = currentJ

        for bj in stride(from: startPos - 1, through: 0, by: -1) {
            let base = query[bj]
            guard base < 4 else { break }

            var numCurr = 0
            var currS: Int64 = -1

            for p in 0..<prev.count {
                let interval = prev[p]
                let ext = BackwardSearch.backwardExt(
                    bwt: bwt,
                    interval: (interval.k, interval.k + interval.s - 1, interval.s),
                    base: Int(base)
                )
                let newS = ext.s

                if newS < minIntv && (interval.n - Int32(bj) + 1) >= minSeedLen {
                    // This interval can't extend further; emit SMEM
                    currentJ = bj
                    let smem = SMEM(
                        k: interval.k,
                        l: interval.k + interval.s - 1,
                        queryBegin: Int32(bj + 1),
                        queryEnd: interval.n + 1
                    )
                    smems.append(smem)
                    break
                }

                if newS >= minIntv && newS != currS {
                    currS = newS
                    if numCurr < prev.count {
                        prev[numCurr] = (ext.k, ext.l, ext.s, Int32(bj), interval.n)
                    }
                    numCurr += 1

                    if p == 0 { break }  // First (longest) interval extended successfully
                }

                // Continue checking shorter intervals
                if p > 0 && newS >= minIntv && newS != currS {
                    currS = newS
                    if numCurr < prev.count {
                        prev[numCurr] = (ext.k, ext.l, ext.s, Int32(bj), interval.n)
                    }
                    numCurr += 1
                }
            }

            if numCurr == 0 { break }
            prev = Array(prev[0..<numCurr])
        }

        // Emit any remaining intervals that weren't emitted during backward extension
        if !prev.isEmpty {
            let interval = prev[0]
            let smemLen = interval.n - interval.m + 1
            if smemLen >= minSeedLen {
                let smem = SMEM(
                    k: interval.k,
                    l: interval.k + interval.s - 1,
                    queryBegin: interval.m,
                    queryEnd: interval.n + 1
                )
                smems.append(smem)
            }
        }

        return (smems, nextPos)
    }
}

// MARK: - BWT helper for count access

extension BWT {
    /// Get cumulative count for base c (accessor for the count tuple).
    @inlinable
    public func count(for c: Int) -> Int64 {
        switch c {
        case 0: return count.0
        case 1: return count.1
        case 2: return count.2
        case 3: return count.3
        case 4: return count.4
        default: return 0
        }
    }

    /// Get cumulative count for base c+1 (next base boundary).
    @inlinable
    public func count(forNext c: Int) -> Int64 {
        switch c {
        case 0: return count.1
        case 1: return count.2
        case 2: return count.3
        case 3: return count.4
        default: return count.4
        }
    }
}
