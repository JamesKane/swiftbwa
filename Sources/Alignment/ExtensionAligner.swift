import BWACore

/// Drives Smith-Waterman extension left and right from seed chains.
/// Reimplements the extension logic from bwa-mem2's `mem_chain2aln_across_reads_V2()`.
public struct ExtensionAligner: Sendable {

    /// Extend a chain into an alignment region using banded Smith-Waterman.
    ///
    /// - Parameters:
    ///   - chain: The seed chain to extend
    ///   - read: The read sequence (2-bit encoded)
    ///   - getReference: Closure to get reference subsequence at (position, length)
    ///   - scoring: Scoring parameters
    /// - Returns: Array of alignment regions (typically one per chain, but may split)
    public static func extend(
        chain: MemChain,
        read: ReadSequence,
        getReference: (Int64, Int) -> [UInt8],
        scoring: ScoringParameters
    ) -> [MemAlnReg] {
        var regions: [MemAlnReg] = []
        let seeds = chain.seeds
        guard !seeds.isEmpty else { return regions }

        let readLen = Int32(read.length)
        let bandWidth = scoring.bandWidth
        let mat = scoring.scoringMatrix()

        // Process each seed as potential extension anchor
        var coveredRegions: [(qb: Int32, qe: Int32, rb: Int64, re: Int64)] = []

        // tmp = max single-event penalty (for subN threshold, matches markSecondaryCore)
        let tmp = max(
            scoring.matchScore + scoring.mismatchPenalty,
            max(scoring.gapOpenPenalty + scoring.gapExtendPenalty,
                scoring.gapOpenPenaltyDeletion + scoring.gapExtendPenaltyDeletion)
        )

        for seedIdx in 0..<seeds.count {
            let seed = seeds[seedIdx]

            // Skip if this seed region is already covered by a previous extension
            let seedQEnd = seed.qbeg + seed.len
            let seedREnd = seed.rbeg + Int64(seed.len)

            var alreadyCovered = false
            var coveringRegIdx = -1
            for (idx, covered) in coveredRegions.enumerated() {
                if seed.qbeg >= covered.qb && seedQEnd <= covered.qe {
                    alreadyCovered = true
                    coveringRegIdx = idx
                    break
                }
            }
            if alreadyCovered {
                // Track sub-optimal score from covered seed (matches bwa-mem2's csub logic).
                // The seed's alignment score estimate is seed.len * matchScore.
                let seedAlnScore = seed.len * scoring.matchScore
                if seedAlnScore > regions[coveringRegIdx].sub {
                    regions[coveringRegIdx].sub = seedAlnScore
                }
                if seedAlnScore >= regions[coveringRegIdx].score - tmp {
                    regions[coveringRegIdx].subN += 1
                }
                continue
            }

            var reg = MemAlnReg()
            reg.w = bandWidth
            reg.rid = chain.rid
            reg.isAlt = chain.isAlt
            reg.seedCov = seed.len
            reg.seedLen0 = Int32(seeds.count)

            // Accumulated h0 chains through extensions (matches bwa-mem2's sc0 pattern).
            // Left extension starts from seed.score; right extension starts from left's
            // local best (or seed.score if left didn't run).
            var accumulatedH0 = seed.score

            // --- Left extension ---
            var leftScore: Int32 = 0
            var leftQLen: Int32 = 0
            var leftTLen: Int32 = 0

            if seed.qbeg > 0 {
                let queryLeft = Array(read.bases[0..<Int(seed.qbeg)].reversed())
                let targetLeftLen = min(Int(seed.qbeg) + Int(bandWidth), Int(seed.rbeg))

                if targetLeftLen > 0 {
                    let targetStart = seed.rbeg - Int64(targetLeftLen)
                    let targetLeft = getReference(targetStart, targetLeftLen).reversed()
                    let targetLeftArr = Array(targetLeft)

                    let result = queryLeft.withUnsafeBufferPointer { qBuf in
                        targetLeftArr.withUnsafeBufferPointer { tBuf in
                            bandedSWExtend(
                                query: qBuf,
                                target: tBuf,
                                scoring: scoring,
                                w: bandWidth,
                                h0: seed.score,
                                scoringMatrix: mat
                            )
                        }
                    }

                    // Update accumulated h0 to left's local best (bwa-mem2: sc0 = a->score)
                    accumulatedH0 = result.score

                    // Clip-vs-extend decision (bwa-mem2 lines 2499-2505)
                    if result.globalScore <= 0
                        || result.globalScore <= result.score - scoring.penClip5 {
                        // CLIP: local alignment is better than paying the clip penalty
                        leftQLen = result.queryEnd
                        leftTLen = result.targetEnd
                        leftScore = result.score
                    } else {
                        // EXTEND TO END: end-to-end alignment preferred
                        leftQLen = seed.qbeg       // extends all the way to position 0
                        leftTLen = result.globalTargetEnd
                        leftScore = result.globalScore
                    }
                }
            }

            // --- Right extension ---
            var rightScore: Int32 = 0
            var rightQLen: Int32 = 0
            var rightTLen: Int32 = 0

            if seedQEnd < readLen {
                let queryRight = Array(read.bases[Int(seedQEnd)..<read.length])
                let targetRightLen = min(Int(readLen - seedQEnd) + Int(bandWidth), 10000)

                if targetRightLen > 0 {
                    let targetRight = getReference(seedREnd, targetRightLen)

                    let result = queryRight.withUnsafeBufferPointer { qBuf in
                        targetRight.withUnsafeBufferPointer { tBuf in
                            bandedSWExtend(
                                query: qBuf,
                                target: tBuf,
                                scoring: scoring,
                                w: bandWidth,
                                h0: accumulatedH0,
                                scoringMatrix: mat
                            )
                        }
                    }

                    // Clip-vs-extend decision (bwa-mem2 lines 2716-2723)
                    if result.globalScore <= 0
                        || result.globalScore <= result.score - scoring.penClip3 {
                        // CLIP: local alignment is better than paying the clip penalty
                        rightQLen = result.queryEnd
                        rightTLen = result.targetEnd
                        rightScore = result.score
                    } else {
                        // EXTEND TO END: end-to-end alignment preferred
                        rightQLen = readLen - seedQEnd  // extends to end of read
                        rightTLen = result.globalTargetEnd
                        rightScore = result.globalScore
                    }
                }
            }

            // Compute final alignment coordinates
            reg.qb = seed.qbeg - leftQLen
            reg.qe = seedQEnd + rightQLen
            reg.rb = seed.rbeg - Int64(leftTLen)
            reg.re = seedREnd + Int64(rightTLen)
            // Score: matches bwa-mem2's truesc accumulation pattern.
            // truesc starts at seed.score, is replaced by leftChosen (which includes
            // h0=seed.score), then adds rightChosen - rightH0 (rightH0 = accumulatedH0).
            var trueScore = Int32(seed.len) * scoring.matchScore
            if leftScore > 0 {
                trueScore = leftScore
            }
            if rightScore > 0 {
                trueScore += rightScore - accumulatedH0
            }
            reg.score = trueScore
            reg.trueScore = trueScore

            // Compute seedCov: total length of seeds fully contained in alignment region
            // Matches bwa-mem2's seedcov computation in mem_chain2aln_across_reads_V2
            var seedCov: Int32 = 0
            for s in seeds {
                if s.qbeg >= reg.qb && s.qbeg + s.len <= reg.qe
                    && s.rbeg >= reg.rb && s.rbeg + Int64(s.len) <= reg.re {
                    seedCov += s.len
                }
            }
            reg.seedCov = seedCov

            coveredRegions.append((reg.qb, reg.qe, reg.rb, reg.re))
            regions.append(reg)
        }

        // Threshold: sub must exceed minSeedLength * matchScore (matches bwa-mem2)
        let subThreshold = scoring.minSeedLength * scoring.matchScore
        for i in 0..<regions.count {
            if regions[i].sub > 0 && regions[i].sub < subThreshold {
                regions[i].sub = 0
            }
        }

        return regions
    }

    /// Tiered SIMD Smith-Waterman extension: try 8-bit SIMD first, fall back to 16-bit.
    /// Matches bwa-mem2's ksw_extend2() dispatch: 8-bit → 16-bit SIMD.
    private static func bandedSWExtend(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        h0: Int32,
        scoringMatrix: [Int8]
    ) -> SWResult {
        if let result = BandedSW8.align(
            query: query, target: target, scoring: scoring, w: w, h0: h0,
            scoringMatrix: scoringMatrix
        ) {
            return result
        }
        return BandedSW16.align(
            query: query, target: target, scoring: scoring, w: w, h0: h0,
            scoringMatrix: scoringMatrix
        )
    }
}
