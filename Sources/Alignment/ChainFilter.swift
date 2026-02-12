import BWACore

/// Filters chains by weight, overlap, and marks primary/secondary.
/// Reimplements `mem_chain_flt()` from bwa-mem2's bwamem.cpp.
public struct ChainFilter: Sendable {

    /// Filter chains: remove low-weight chains and those significantly overlapping heavier chains.
    public static func filter(
        chains: inout [MemChain],
        scoring: ScoringParameters
    ) {
        guard !chains.isEmpty else { return }

        let dropRatio = scoring.chainDropRatio
        let minWeight = scoring.minChainWeight > 0 ? scoring.minChainWeight : scoring.minSeedLength

        // Remove chains below minimum weight
        chains.removeAll { $0.weight < minWeight }
        guard !chains.isEmpty else { return }

        // Sort by weight descending
        chains.sort { $0.weight > $1.weight }

        // Mark kept/removed based on overlap with heavier chains
        for i in 0..<chains.count {
            if chains[i].kept == 0 { continue }

            let iBegin = chains[i].seeds.first?.qbeg ?? 0
            let iEnd = chains[i].seeds.last.map { $0.qbeg + $0.len } ?? 0

            for j in (i + 1)..<chains.count {
                if chains[j].kept == 0 { continue }

                let jBegin = chains[j].seeds.first?.qbeg ?? 0
                let jEnd = chains[j].seeds.last.map { $0.qbeg + $0.len } ?? 0

                // Calculate overlap
                let overlapBegin = max(iBegin, jBegin)
                let overlapEnd = min(iEnd, jEnd)

                // Don't suppress primary chain (j) when overlapping ALT chain (i)
                if overlapBegin < overlapEnd
                    && (!chains[i].isAlt || chains[j].isAlt) {
                    let overlap = overlapEnd - overlapBegin
                    let jLen = jEnd - jBegin

                    if jLen > 0 && Float(overlap) / Float(jLen) > dropRatio &&
                       Float(chains[j].weight) < Float(chains[i].weight) * dropRatio {
                        chains[j].kept = 0
                    }
                }
            }
        }

        // Remove unkept chains
        chains.removeAll { $0.kept == 0 }
    }

    /// Mark secondary alignment regions based on overlap with primary regions.
    public static func markSecondary(
        regions: inout [MemAlnReg],
        maskLevel: Float,
        readId: UInt64 = 0
    ) {
        guard regions.count > 1 else { return }

        // Assign hashes for deterministic tiebreaking (matches bwa-mem2's mem_mark_primary_se)
        for i in 0..<regions.count {
            regions[i].hash = hash64(readId &+ UInt64(i))
        }

        // Sort by score descending, then hash for deterministic tiebreaking
        regions.sort { a, b in
            if a.score != b.score { return a.score > b.score }
            return a.hash < b.hash
        }

        for i in 0..<regions.count {
            if regions[i].secondary >= 0 { continue }

            for j in (i + 1)..<regions.count {
                if regions[j].secondary >= 0 { continue }

                // Check overlap on query
                let qOverlapBegin = max(regions[i].qb, regions[j].qb)
                let qOverlapEnd = min(regions[i].qe, regions[j].qe)

                if qOverlapBegin < qOverlapEnd {
                    let overlap = Float(qOverlapEnd - qOverlapBegin)
                    let minLen = Float(min(regions[i].qe - regions[i].qb, regions[j].qe - regions[j].qb))

                    if minLen > 0 && overlap / minLen > maskLevel {
                        regions[j].secondary = Int32(i)
                    }
                }
            }
        }
    }

    // MARK: - ALT-Aware Two-Phase Secondary Marking

    /// Two-phase secondary marking for ALT-aware alignment.
    /// Implements bwa-mem2's `mem_mark_primary_se()`:
    /// - Phase 1: Mark among ALL hits (ALT + primary), recording `secondaryAll`
    /// - Phase 2: Re-mark among primary-only subset so ALT hits cannot suppress primaries
    public static func markSecondaryALT(
        regions: inout [MemAlnReg],
        maskLevel: Float,
        scoring: ScoringParameters,
        readId: UInt64 = 0
    ) {
        guard regions.count > 1 else { return }
        let n = regions.count

        // Initialize (matches bwa-mem2's mem_mark_primary_se initialization)
        for i in 0..<n {
            regions[i].sub = 0
            regions[i].altSc = 0
            regions[i].secondary = -1
            regions[i].secondaryAll = -1
            regions[i].hash = hash64(readId &+ UInt64(i))
        }

        let nPri = regions.filter { !$0.isAlt }.count

        // Phase 1: Sort ALL hits by score desc, then isAlt asc (primary before ALT)
        regions.sort { a, b in
            if a.score != b.score { return a.score > b.score }
            if a.isAlt != b.isAlt { return !a.isAlt }
            return a.hash < b.hash
        }
        markSecondaryCore(regions: &regions, count: n, maskLevel: maskLevel, scoring: scoring)

        // Record Phase 1 rankings
        for i in 0..<n {
            regions[i].secondaryAll = Int32(i)
            // Track altSc: if primary hit is secondary to an ALT hit
            if !regions[i].isAlt && regions[i].secondary >= 0
                && regions[Int(regions[i].secondary)].isAlt {
                regions[i].altSc = regions[Int(regions[i].secondary)].score
            }
        }

        // Phase 2: If mixed ALT+primary, re-rank among primary-only
        if nPri > 0 && nPri < n {
            // Re-sort: primary first, then by score
            regions.sort { a, b in
                if a.isAlt != b.isAlt { return !a.isAlt }
                if a.score != b.score { return a.score > b.score }
                return a.hash < b.hash
            }

            // Build re-index mapping: old secondaryAll -> new position
            var reindex = [Int](repeating: 0, count: n)
            for i in 0..<n {
                reindex[Int(regions[i].secondaryAll)] = i
            }

            // Update secondary pointers and mark ALT secondaries
            for i in 0..<n {
                if regions[i].secondary >= 0 {
                    regions[i].secondaryAll = Int32(reindex[Int(regions[i].secondary)])
                    if regions[i].isAlt {
                        regions[i].secondary = Int32.max
                    }
                } else {
                    regions[i].secondaryAll = -1
                }
            }

            // Re-mark primary-only subset
            for i in 0..<nPri {
                regions[i].sub = 0
                regions[i].secondary = -1
            }
            markSecondaryCore(
                regions: &regions, count: nPri, maskLevel: maskLevel, scoring: scoring
            )
        } else {
            // No mixed: secondaryAll = secondary
            for i in 0..<n {
                regions[i].secondaryAll = regions[i].secondary
            }
        }
    }

    /// Core overlap-based secondary marking with ALT-aware sub_n tracking.
    /// Operates on the first `count` elements of the regions array.
    private static func markSecondaryCore(
        regions: inout [MemAlnReg],
        count: Int,
        maskLevel: Float,
        scoring: ScoringParameters
    ) {
        guard count > 1 else { return }

        // tmp = max single-event penalty (for subN threshold)
        let tmp = max(
            scoring.matchScore + scoring.mismatchPenalty,
            max(scoring.gapOpenPenalty + scoring.gapExtendPenalty,
                scoring.gapOpenPenaltyDeletion + scoring.gapExtendPenaltyDeletion)
        )

        var primaries: [Int] = [0]
        for i in 1..<count {
            var found = false
            for pi in primaries {
                let overlapBegin = max(regions[pi].qb, regions[i].qb)
                let overlapEnd = min(regions[pi].qe, regions[i].qe)
                if overlapBegin < overlapEnd {
                    let overlap = Float(overlapEnd - overlapBegin)
                    let minLen = Float(min(
                        regions[pi].qe - regions[pi].qb,
                        regions[i].qe - regions[i].qb
                    ))
                    if minLen > 0 && overlap / minLen > maskLevel {
                        if regions[pi].sub == 0 { regions[pi].sub = regions[i].score }
                        if regions[pi].score - regions[i].score <= tmp
                            && (regions[pi].isAlt || !regions[i].isAlt) {
                            regions[pi].subN += 1
                        }
                        regions[i].secondary = Int32(pi)
                        found = true
                        break
                    }
                }
            }
            if !found {
                primaries.append(i)
            }
        }
    }
}
