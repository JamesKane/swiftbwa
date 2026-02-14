import BWACore

/// Filters chains by weight, overlap, and marks primary/secondary.
/// Reimplements `mem_chain_flt()` from bwa-mem2's bwamem.cpp.
public struct ChainFilter: Sendable {

    /// Filter chains: remove low-weight chains and those significantly overlapping heavier chains.
    /// Matches bwa-mem2's `mem_chain_flt()` (bwamem.cpp:506-624), including:
    /// - Symmetric overlap formula: `overlap >= min(len_i, len_j) * mask_level`
    /// - Absolute weight difference check: `weight_diff >= minSeedLen * 2`
    /// - First-shadowed-chain recovery: one suppressed chain per kept chain is
    ///   recovered (kept=1) so it produces a region for accurate MAPQ/sub computation.
    /// Compute chain weight as min(non-overlapping query coverage, non-overlapping reference coverage).
    /// Matches bwa-mem2's `mem_chain_weight()` (bwamem.cpp:429-448).
    /// Seeds must already be sorted by reference position.
    private static func chainWeight(_ chain: MemChain) -> Int32 {
        // Pass 1: non-overlapping query coverage
        var w: Int32 = 0
        var end: Int64 = 0
        for seed in chain.seeds {
            let sb = Int64(seed.qbeg)
            let se = sb + Int64(seed.len)
            if sb >= end {
                w += seed.len
            } else if se > end {
                w += Int32(se - end)
            }
            if se > end { end = se }
        }
        let tmp = w

        // Pass 2: non-overlapping reference coverage
        w = 0
        end = 0
        for seed in chain.seeds {
            let sb = seed.rbeg
            let se = sb + Int64(seed.len)
            if sb >= end {
                w += seed.len
            } else if se > end {
                w += Int32(se - end)
            }
            if se > end { end = se }
        }

        w = min(w, tmp)
        return w < (1 << 30) ? w : (1 << 30) - 1
    }

    public static func filter(
        chains: inout [MemChain],
        scoring: ScoringParameters
    ) {
        guard !chains.isEmpty else { return }

        let dropRatio = scoring.chainDropRatio
        let maskLevel = scoring.maskLevel
        let minWeight = scoring.minChainWeight > 0 ? scoring.minChainWeight : scoring.minSeedLength

        // Recompute chain weight properly (bwa-mem2 line 515: c->w = mem_chain_weight(c))
        for i in 0..<chains.count {
            chains[i].weight = chainWeight(chains[i])
        }

        // Remove chains below minimum weight
        chains.removeAll { $0.weight < minWeight }
        guard !chains.isEmpty else { return }

        // Sort by weight descending
        chains.sort { $0.weight > $1.weight }

        // Initialize kept/first for all chains
        for i in 0..<chains.count {
            chains[i].kept = 0
            chains[i].first = -1
        }

        // Heaviest chain is always kept
        chains[0].kept = 3
        var keptIndices: [Int] = [0]

        // Process chains from second-heaviest to lightest (bwa-mem2 lines 565-589)
        for i in 1..<chains.count {
            let iBegin = chains[i].seeds.first?.qbeg ?? 0
            let iEnd = chains[i].seeds.last.map { $0.qbeg + $0.len } ?? 0
            let iLen = iEnd - iBegin

            var largeOvlp = false
            var suppressed = false

            for ki in keptIndices {
                let kBegin = chains[ki].seeds.first?.qbeg ?? 0
                let kEnd = chains[ki].seeds.last.map { $0.qbeg + $0.len } ?? 0
                let kLen = kEnd - kBegin

                let overlapBegin = max(iBegin, kBegin)
                let overlapEnd = min(iEnd, kEnd)

                // Don't let ALT kept chain suppress primary candidate
                if overlapBegin < overlapEnd
                    && (!chains[ki].isAlt || chains[i].isAlt)
                {
                    let overlap = overlapEnd - overlapBegin
                    let minLen = min(iLen, kLen)

                    // Significant overlap: overlap >= min_l * mask_level
                    // Use float comparison to match C's int * float promotion
                    if minLen > 0
                        && Float(overlap) >= Float(minLen) * maskLevel
                        && minLen < scoring.maxChainGap
                    {
                        largeOvlp = true
                        // Record first shadowed chain for this kept chain
                        if chains[ki].first < 0 {
                            chains[ki].first = Int32(i)
                        }
                        // Suppress if weight ratio AND absolute difference both exceeded
                        if Float(chains[i].weight) < Float(chains[ki].weight) * dropRatio
                            && chains[ki].weight - chains[i].weight >= scoring.minSeedLength << 1
                        {
                            suppressed = true
                            break
                        }
                    }
                }
            }

            if !suppressed {
                keptIndices.append(i)
                chains[i].kept = largeOvlp ? 2 : 3
            }
            // else chains[i].kept stays 0 (will be removed)
        }

        // Recover first-shadowed chains for MAPQ (bwa-mem2 lines 591-595)
        for ki in keptIndices {
            let firstIdx = Int(chains[ki].first)
            if firstIdx >= 0 && firstIdx < chains.count {
                chains[firstIdx].kept = 1
            }
        }

        // Remove unkept chains (max_chain_extend is effectively unlimited)
        chains.removeAll { $0.kept == 0 }
    }

    /// Mark secondary alignment regions based on overlap with primary regions.
    /// Delegates to `markSecondaryCore` so that `sub` and `subN` are populated
    /// (critical for correct MAPQ computation).
    public static func markSecondary(
        regions: inout [MemAlnReg],
        maskLevel: Float,
        scoring: ScoringParameters = ScoringParameters(),
        readId: UInt64 = 0
    ) {
        guard regions.count > 1 else { return }

        for i in 0..<regions.count {
            regions[i].sub = 0
            regions[i].subN = 0
            regions[i].secondary = -1
            regions[i].hash = hash64(readId &+ UInt64(i))
        }

        // Sort by score descending, then hash for deterministic tiebreaking
        regions.sort { a, b in
            if a.score != b.score { return a.score > b.score }
            return a.hash < b.hash
        }

        markSecondaryCore(
            regions: &regions, count: regions.count,
            maskLevel: maskLevel, scoring: scoring
        )
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
