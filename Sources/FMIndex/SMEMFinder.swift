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
    /// Find all SMEMs with a custom minimum SA interval threshold.
    ///
    /// Used for re-seeding: `minIntv > 1` causes the search to ignore seeds with
    /// fewer occurrences, finding additional shorter but more specific seeds.
    public static func findAllSMEMs(
        query: [UInt8],
        bwt: BWT,
        minSeedLen: Int32,
        minIntv: Int64 = 1
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
                minIntv: minIntv
            )
            allSMEMs.append(contentsOf: smems)
            pos = nextPos
        }

        allSMEMs.sort {
            if $0.queryBegin != $1.queryBegin {
                return $0.queryBegin < $1.queryBegin
            }
            return $0.length > $1.length
        }

        return allSMEMs
    }

    /// Bidirectional SMEM interval used during search.
    /// Tracks both forward (k) and reverse (l) SA intervals plus query span.
    public struct BidiInterval {
        public var k: Int64    // Forward SA interval start
        public var l: Int64    // Reverse SA interval position
        public var s: Int64    // Interval size (same in both directions)
        public var m: Int32    // Query begin (inclusive)
        public var n: Int32    // Query end (inclusive, i.e., last matched position)
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

        // Initialize bidirectional SA interval for the first character
        // Forward: k = count[a], Reverse: l = count[3-a], s = count[a+1] - count[a]
        var smem = BidiInterval(
            k: bwt.count(for: Int(a)),
            l: bwt.count(for: Int(3 - a)),
            s: bwt.count(forNext: Int(a)) - bwt.count(for: Int(a)),
            m: Int32(startPos),
            n: Int32(startPos)
        )

        // Track previous intervals (candidates for backward extension)
        var prevArray: [BidiInterval] = []
        var numPrev = 0

        // === Forward phase: extend rightward (lines 537-575) ===
        var j = startPos + 1
        while j < readLength {
            let base = query[j]
            nextPos = j + 1

            guard base < 4 else { break }

            // Forward extension = backward extension with swapped k/l and complement base
            let swapped = BidiInterval(k: smem.l, l: smem.k, s: smem.s, m: smem.m, n: smem.n)
            let ext = BackwardSearch.backwardExt(
                bwt: bwt,
                interval: (swapped.k, swapped.l, swapped.s),
                base: Int(3 - base)
            )
            // Swap result back
            let newSmem = BidiInterval(k: ext.l, l: ext.k, s: ext.s, m: smem.m, n: Int32(j))

            // If interval size changed, record previous as candidate
            // (using bitmask trick from C++: prevArray[numPrev] = smem; numPrev += (newS != oldS))
            if newSmem.s != smem.s {
                prevArray.append(smem)
                numPrev += 1
            }

            // If new interval is too small, stop forward extension
            if newSmem.s < minIntv {
                nextPos = j
                break
            }

            smem = newSmem
            j += 1
        }

        // Record the final forward interval if it meets threshold
        if smem.s >= minIntv {
            prevArray.append(smem)
            numPrev += 1
        }

        // Reverse prevArray so longest match is first (lines 587-592)
        prevArray.reverse()

        // === Backward phase: extend leftward from each candidate (lines 596-655) ===
        for bj in stride(from: startPos - 1, through: 0, by: -1) {
            let base = query[bj]
            guard base < 4 else { break }

            var numCurr = 0
            var currS: Int64 = -1

            // First loop: p = 0..<numPrev (lines 607-629)
            var p = 0
            while p < numPrev {
                let interval = prevArray[p]
                let ext = BackwardSearch.backwardExt(
                    bwt: bwt,
                    interval: (interval.k, interval.l, interval.s),
                    base: Int(base)
                )
                let newInterval = BidiInterval(k: ext.k, l: ext.l, s: ext.s, m: Int32(bj), n: interval.n)

                // CONDITION 1: Can't extend further AND meets min length
                if newInterval.s < minIntv && (interval.n - interval.m + 1) >= minSeedLen {
                    // Emit the PREVIOUS interval (not the new one)
                    let emitSmem = SMEM(
                        k: interval.k,
                        l: interval.k + interval.s - 1,
                        queryBegin: interval.m,
                        queryEnd: interval.n + 1
                    )
                    smems.append(emitSmem)
                    break
                }

                // CONDITION 2: Can extend AND hasn't seen this size yet
                if newInterval.s >= minIntv && newInterval.s != currS {
                    currS = newInterval.s
                    if numCurr < prevArray.count {
                        prevArray[numCurr] = newInterval
                    } else {
                        prevArray.append(newInterval)
                    }
                    numCurr += 1
                    break  // Break first loop (line 628)
                }
                p += 1
            }

            // Second loop: p+1..<numPrev (lines 631-649)
            p += 1
            while p < numPrev {
                let interval = prevArray[p]
                let ext = BackwardSearch.backwardExt(
                    bwt: bwt,
                    interval: (interval.k, interval.l, interval.s),
                    base: Int(base)
                )
                let newInterval = BidiInterval(k: ext.k, l: ext.l, s: ext.s, m: Int32(bj), n: interval.n)

                if newInterval.s >= minIntv && newInterval.s != currS {
                    currS = newInterval.s
                    if numCurr < prevArray.count {
                        prevArray[numCurr] = newInterval
                    } else {
                        prevArray.append(newInterval)
                    }
                    numCurr += 1
                }
                p += 1
            }

            numPrev = numCurr
            if numCurr == 0 {
                break
            }
        }

        // Emit any remaining intervals (lines 656-664)
        if numPrev != 0 {
            let interval = prevArray[0]
            if (interval.n - interval.m + 1) >= minSeedLen {
                let emitSmem = SMEM(
                    k: interval.k,
                    l: interval.k + interval.s - 1,
                    queryBegin: interval.m,
                    queryEnd: interval.n + 1
                )
                smems.append(emitSmem)
            }
        }

        return (smems, nextPos)
    }

    // MARK: - Arena-backed overloads

    /// Arena-backed overload: appends to caller's buffers, zero per-call allocations.
    public static func findSMEMsAtPosition(
        query: [UInt8],
        bwt: BWT,
        startPos: Int,
        minSeedLen: Int32,
        minIntv: Int64,
        smemOutput: inout ArenaBuffer<SMEM>,
        prevBuffer: inout ArenaBuffer<BidiInterval>
    ) -> Int {
        let readLength = query.count
        var nextPos = startPos + 1

        let a = query[startPos]
        guard a < 4 else {
            return nextPos
        }

        var smem = BidiInterval(
            k: bwt.count(for: Int(a)),
            l: bwt.count(for: Int(3 - a)),
            s: bwt.count(forNext: Int(a)) - bwt.count(for: Int(a)),
            m: Int32(startPos),
            n: Int32(startPos)
        )

        var numPrev = 0

        // === Forward phase: extend rightward ===
        var j = startPos + 1
        while j < readLength {
            let base = query[j]
            nextPos = j + 1

            guard base < 4 else { break }

            let swapped = BidiInterval(k: smem.l, l: smem.k, s: smem.s, m: smem.m, n: smem.n)
            let ext = BackwardSearch.backwardExt(
                bwt: bwt,
                interval: (swapped.k, swapped.l, swapped.s),
                base: Int(3 - base)
            )
            let newSmem = BidiInterval(k: ext.l, l: ext.k, s: ext.s, m: smem.m, n: Int32(j))

            if newSmem.s != smem.s {
                prevBuffer.append(smem)
                numPrev += 1
            }

            if newSmem.s < minIntv {
                nextPos = j
                break
            }

            smem = newSmem
            j += 1
        }

        if smem.s >= minIntv {
            prevBuffer.append(smem)
            numPrev += 1
        }

        // Reverse prevBuffer portion in-place
        let base0 = prevBuffer.count - numPrev
        var lo = base0, hi = prevBuffer.count - 1
        while lo < hi {
            let tmp = prevBuffer[lo]
            prevBuffer[lo] = prevBuffer[hi]
            prevBuffer[hi] = tmp
            lo += 1; hi -= 1
        }

        // === Backward phase: extend leftward from each candidate ===
        for bj in stride(from: startPos - 1, through: 0, by: -1) {
            let base = query[bj]
            guard base < 4 else { break }

            var numCurr = 0
            var currS: Int64 = -1

            var p = 0
            while p < numPrev {
                let interval = prevBuffer[base0 + p]
                let ext = BackwardSearch.backwardExt(
                    bwt: bwt,
                    interval: (interval.k, interval.l, interval.s),
                    base: Int(base)
                )
                let newInterval = BidiInterval(k: ext.k, l: ext.l, s: ext.s, m: Int32(bj), n: interval.n)

                if newInterval.s < minIntv && (interval.n - interval.m + 1) >= minSeedLen {
                    let emitSmem = SMEM(
                        k: interval.k,
                        l: interval.k + interval.s - 1,
                        queryBegin: interval.m,
                        queryEnd: interval.n + 1
                    )
                    smemOutput.append(emitSmem)
                    break
                }

                if newInterval.s >= minIntv && newInterval.s != currS {
                    currS = newInterval.s
                    prevBuffer[base0 + numCurr] = newInterval
                    numCurr += 1
                    break
                }
                p += 1
            }

            p += 1
            while p < numPrev {
                let interval = prevBuffer[base0 + p]
                let ext = BackwardSearch.backwardExt(
                    bwt: bwt,
                    interval: (interval.k, interval.l, interval.s),
                    base: Int(base)
                )
                let newInterval = BidiInterval(k: ext.k, l: ext.l, s: ext.s, m: Int32(bj), n: interval.n)

                if newInterval.s >= minIntv && newInterval.s != currS {
                    currS = newInterval.s
                    prevBuffer[base0 + numCurr] = newInterval
                    numCurr += 1
                }
                p += 1
            }

            numPrev = numCurr
            if numCurr == 0 { break }
        }

        if numPrev != 0 {
            let interval = prevBuffer[base0]
            if (interval.n - interval.m + 1) >= minSeedLen {
                let emitSmem = SMEM(
                    k: interval.k,
                    l: interval.k + interval.s - 1,
                    queryBegin: interval.m,
                    queryEnd: interval.n + 1
                )
                smemOutput.append(emitSmem)
            }
        }

        return nextPos
    }

    /// Arena-backed: allocates SMEM + prev buffers from arena, returns ArenaBuffer<SMEM>.
    public static func findAllSMEMs(
        query: [UInt8],
        bwt: BWT,
        minSeedLen: Int32,
        minIntv: Int64 = 1,
        arena: inout ReadArena
    ) -> ArenaBuffer<SMEM> {
        let readLen = query.count
        guard readLen > 0 else {
            return ArenaBuffer(base: arena.allocate(SMEM.self, count: 1), capacity: 1)
        }

        let smemCapacity = max(readLen * 16, 4096)
        var smems = ArenaBuffer<SMEM>(
            base: arena.allocate(SMEM.self, count: smemCapacity),
            capacity: smemCapacity
        )
        var prevBuf = ArenaBuffer<BidiInterval>(
            base: arena.allocate(BidiInterval.self, count: readLen),
            capacity: readLen
        )

        var pos = 0
        while pos < readLen {
            prevBuf.removeAll()
            pos = findSMEMsAtPosition(
                query: query, bwt: bwt, startPos: pos,
                minSeedLen: minSeedLen, minIntv: minIntv,
                smemOutput: &smems, prevBuffer: &prevBuf
            )
        }

        smems.sort {
            if $0.queryBegin != $1.queryBegin { return $0.queryBegin < $1.queryBegin }
            return $0.length > $1.length
        }
        return smems
    }
}
