import BWACore

/// Pre-computed extension plan for a single seed, used by GPU batch path.
/// Contains pre-fetched, reversed sequences for left and forward sequences for right.
public struct SeedExtensionPlan: Sendable {
    public let chainIndex: Int
    public let seedIndexInChain: Int
    public let seed: MemSeed
    public let seedQEnd: Int32
    public let seedREnd: Int64
    public let readLen: Int32
    public let chainRid: Int32
    public let chainIsAlt: Bool
    public let chainFracRep: Float
    public let chainSeedCount: Int
    /// All seeds in the parent chain (for seedCov computation).
    public let chainSeeds: [MemSeed]
    /// Reversed query bases for left extension (empty if seed.qbeg == 0).
    public let leftQuery: [UInt8]
    /// Reversed reference bases for left extension (empty if no left extension needed).
    public let leftTarget: [UInt8]
    /// Forward query bases for right extension (empty if seedQEnd == readLen).
    public let rightQuery: [UInt8]
    /// Forward reference bases for right extension (empty if no right extension needed).
    public let rightTarget: [UInt8]
}

/// Drives Smith-Waterman extension left and right from seed chains.
/// Reimplements the extension logic from bwa-mem2's `mem_chain2aln_across_reads_V2()`.
public struct ExtensionAligner: Sendable {

    /// SIMD4-vectorized seedCov: total length of seeds fully contained in a region.
    /// Query check uses SIMD4<Int32> (4 seeds per NEON instruction);
    /// reference check remains scalar (Int64 needs SIMD2, limited benefit).
    @inline(__always)
    private static func computeSeedCov(
        seeds: SeedSoA,
        qb: Int32, qe: Int32, rb: Int64, re: Int64
    ) -> Int32 {
        let sq = seeds.qbegs, sl = seeds.lens, sr = seeds.rbegs
        let n = seeds.count
        var seedCov: Int32 = 0

        let regQb = SIMD4<Int32>(repeating: qb)
        let regQe = SIMD4<Int32>(repeating: qe)
        var i = 0
        while i + 4 <= n {
            let sqV = SIMD4(sq[i], sq[i+1], sq[i+2], sq[i+3])
            let slV = SIMD4(sl[i], sl[i+1], sl[i+2], sl[i+3])
            let sqEnd = sqV &+ slV
            let qMask = (sqV .>= regQb) .& (sqEnd .<= regQe)
            // Reference check requires Int64 — scalar for masked lanes
            for j in 0..<4 where qMask[j] {
                let idx = i + j
                if sr[idx] >= rb && sr[idx] + Int64(sl[idx]) <= re {
                    seedCov += sl[idx]
                }
            }
            i += 4
        }
        // Scalar tail
        while i < n {
            if sq[i] >= qb && sq[i] + sl[i] <= qe
                && sr[i] >= rb && sr[i] + Int64(sl[i]) <= re {
                seedCov += sl[i]
            }
            i += 1
        }
        return seedCov
    }

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
        guard !chain.seeds.isEmpty else { return regions }

        let soa = SeedSoA(from: chain.seeds)
        defer { soa.deallocate() }

        let readLen = Int32(read.length)
        let bandWidth = scoring.bandWidth
        let mat = scoring.scoringMatrix()

        // tmp = max single-event penalty (for subN threshold, matches markSecondaryCore)
        let tmp = max(
            scoring.matchScore + scoring.mismatchPenalty,
            max(scoring.gapOpenPenalty + scoring.gapExtendPenalty,
                scoring.gapOpenPenaltyDeletion + scoring.gapExtendPenaltyDeletion)
        )

        // Within-chain coverage: track query ranges of already-extended regions
        var coveredRegions: [(qb: Int32, qe: Int32)] = []

        for seedIdx in 0..<soa.count {
            let seedQbeg = soa.qbegs[seedIdx]
            let seedLen = soa.lens[seedIdx]
            let seedRbeg = soa.rbegs[seedIdx]
            let seedScore = soa.scores[seedIdx]

            let seedQEnd = seedQbeg + seedLen
            let seedREnd = seedRbeg + Int64(seedLen)

            // Within-chain coverage check: simple query containment
            var coveringRegIdx = -1
            for (idx, cov) in coveredRegions.enumerated() {
                if seedQbeg >= cov.qb && seedQEnd <= cov.qe {
                    coveringRegIdx = idx
                    break
                }
            }

            if coveringRegIdx >= 0 {
                // Seed is covered within this chain — update sub/subN
                let seedAlnScore = seedLen * scoring.matchScore
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
            reg.fracRep = chain.fracRep
            reg.seedCov = seedLen
            reg.seedLen0 = seedLen

            // Accumulated h0 chains through extensions (matches bwa-mem2's sc0 pattern).
            // Left extension starts from seed.score; right extension starts from left's
            // local best (or seed.score if left didn't run).
            var accumulatedH0 = seedScore

            // --- Left extension ---
            var leftScore: Int32 = 0
            var leftQLen: Int32 = 0
            var leftTLen: Int32 = 0

            if seedQbeg > 0 {
                let queryLeft = Array(read.bases[0..<Int(seedQbeg)].reversed())
                let targetLeftLen = min(Int(seedQbeg) + Int(bandWidth), Int(seedRbeg))

                if targetLeftLen > 0 {
                    let targetStart = seedRbeg - Int64(targetLeftLen)
                    let targetLeft = getReference(targetStart, targetLeftLen).reversed()
                    let targetLeftArr = Array(targetLeft)

                    let result = queryLeft.withUnsafeBufferPointer { qBuf in
                        targetLeftArr.withUnsafeBufferPointer { tBuf in
                            bandedSWExtend(
                                query: qBuf,
                                target: tBuf,
                                scoring: scoring,
                                w: bandWidth,
                                h0: seedScore,
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
                        leftQLen = seedQbeg       // extends all the way to position 0
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
            reg.qb = seedQbeg - leftQLen
            reg.qe = seedQEnd + rightQLen
            reg.rb = seedRbeg - Int64(leftTLen)
            reg.re = seedREnd + Int64(rightTLen)
            // Score: matches bwa-mem2's truesc accumulation pattern.
            // truesc starts at seed.score, is replaced by leftChosen (which includes
            // h0=seed.score), then adds rightChosen - rightH0 (rightH0 = accumulatedH0).
            var trueScore = Int32(seedLen) * scoring.matchScore
            if leftScore > 0 {
                trueScore = leftScore
            }
            if rightScore > 0 {
                trueScore += rightScore - accumulatedH0
            }
            reg.score = trueScore
            reg.trueScore = trueScore

            reg.seedCov = computeSeedCov(seeds: soa, qb: reg.qb, qe: reg.qe, rb: reg.rb, re: reg.re)

            coveredRegions.append((qb: reg.qb, qe: reg.qe))
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

    /// Extend a chain using pre-allocated workspace buffers to avoid per-seed heap allocations.
    ///
    /// - Parameters:
    ///   - chain: The seed chain to extend
    ///   - read: The read sequence (2-bit encoded)
    ///   - getReference: Closure that writes reference bases into a buffer and returns count written
    ///   - scoring: Scoring parameters
    ///   - scoringMatrix: Pre-computed 25-element scoring matrix
    ///   - queryBuf: Pre-allocated buffer of at least `read.length` bytes
    ///   - targetBuf: Pre-allocated buffer of at least 10200 bytes
    public static func extend(
        chain: MemChain,
        read: ReadSequence,
        getReference: (Int64, Int, UnsafeMutablePointer<UInt8>) -> Int,
        scoring: ScoringParameters,
        scoringMatrix: [Int8],
        queryBuf: UnsafeMutablePointer<UInt8>,
        targetBuf: UnsafeMutablePointer<UInt8>
    ) -> [MemAlnReg] {
        var regions: [MemAlnReg] = []
        guard !chain.seeds.isEmpty else { return regions }

        let soa = SeedSoA(from: chain.seeds)
        defer { soa.deallocate() }

        let readLen = Int32(read.length)
        let bandWidth = scoring.bandWidth
        let mat = scoringMatrix

        let tmp = max(
            scoring.matchScore + scoring.mismatchPenalty,
            max(scoring.gapOpenPenalty + scoring.gapExtendPenalty,
                scoring.gapOpenPenaltyDeletion + scoring.gapExtendPenaltyDeletion)
        )

        var coveredRegions: [(qb: Int32, qe: Int32)] = []

        for seedIdx in 0..<soa.count {
            let seedQbeg = soa.qbegs[seedIdx]
            let seedLen = soa.lens[seedIdx]
            let seedRbeg = soa.rbegs[seedIdx]
            let seedScore = soa.scores[seedIdx]

            let seedQEnd = seedQbeg + seedLen
            let seedREnd = seedRbeg + Int64(seedLen)

            var coveringRegIdx = -1
            for (idx, cov) in coveredRegions.enumerated() {
                if seedQbeg >= cov.qb && seedQEnd <= cov.qe {
                    coveringRegIdx = idx
                    break
                }
            }

            if coveringRegIdx >= 0 {
                let seedAlnScore = seedLen * scoring.matchScore
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
            reg.fracRep = chain.fracRep
            reg.seedCov = seedLen
            reg.seedLen0 = seedLen

            var accumulatedH0 = seedScore

            // --- Left extension ---
            var leftScore: Int32 = 0
            var leftQLen: Int32 = 0
            var leftTLen: Int32 = 0

            if seedQbeg > 0 {
                let queryLeftLen = Int(seedQbeg)
                // Write reversed query into queryBuf
                for i in 0..<queryLeftLen {
                    queryBuf[i] = read.bases[queryLeftLen - 1 - i]
                }

                let targetLeftLen = min(queryLeftLen + Int(bandWidth), Int(seedRbeg))
                if targetLeftLen > 0 {
                    let targetStart = seedRbeg - Int64(targetLeftLen)
                    let actualTargetLen = getReference(targetStart, targetLeftLen, targetBuf)

                    // Reverse targetBuf in-place
                    var lo = 0; var hi = actualTargetLen - 1
                    while lo < hi {
                        let t = targetBuf[lo]; targetBuf[lo] = targetBuf[hi]; targetBuf[hi] = t
                        lo += 1; hi -= 1
                    }

                    let qBuf = UnsafeBufferPointer(start: queryBuf, count: queryLeftLen)
                    let tBuf = UnsafeBufferPointer(start: targetBuf, count: actualTargetLen)
                    let result = bandedSWExtend(
                        query: qBuf, target: tBuf,
                        scoring: scoring, w: bandWidth,
                        h0: seedScore, scoringMatrix: mat
                    )

                    accumulatedH0 = result.score

                    if result.globalScore <= 0
                        || result.globalScore <= result.score - scoring.penClip5 {
                        leftQLen = result.queryEnd
                        leftTLen = result.targetEnd
                        leftScore = result.score
                    } else {
                        leftQLen = seedQbeg
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
                let queryRightLen = read.length - Int(seedQEnd)
                // Write forward query into queryBuf
                for i in 0..<queryRightLen {
                    queryBuf[i] = read.bases[Int(seedQEnd) + i]
                }

                let targetRightLen = min(Int(readLen - seedQEnd) + Int(bandWidth), 10000)
                if targetRightLen > 0 {
                    let actualTargetLen = getReference(seedREnd, targetRightLen, targetBuf)

                    let qBuf = UnsafeBufferPointer(start: queryBuf, count: queryRightLen)
                    let tBuf = UnsafeBufferPointer(start: targetBuf, count: actualTargetLen)
                    let result = bandedSWExtend(
                        query: qBuf, target: tBuf,
                        scoring: scoring, w: bandWidth,
                        h0: accumulatedH0, scoringMatrix: mat
                    )

                    if result.globalScore <= 0
                        || result.globalScore <= result.score - scoring.penClip3 {
                        rightQLen = result.queryEnd
                        rightTLen = result.targetEnd
                        rightScore = result.score
                    } else {
                        rightQLen = readLen - seedQEnd
                        rightTLen = result.globalTargetEnd
                        rightScore = result.globalScore
                    }
                }
            }

            reg.qb = seedQbeg - leftQLen
            reg.qe = seedQEnd + rightQLen
            reg.rb = seedRbeg - Int64(leftTLen)
            reg.re = seedREnd + Int64(rightTLen)

            var trueScore = Int32(seedLen) * scoring.matchScore
            if leftScore > 0 {
                trueScore = leftScore
            }
            if rightScore > 0 {
                trueScore += rightScore - accumulatedH0
            }
            reg.score = trueScore
            reg.trueScore = trueScore

            reg.seedCov = computeSeedCov(seeds: soa, qb: reg.qb, qe: reg.qe, rb: reg.rb, re: reg.re)

            coveredRegions.append((qb: reg.qb, qe: reg.qe))
            regions.append(reg)
        }

        let subThreshold = scoring.minSeedLength * scoring.matchScore
        for i in 0..<regions.count {
            if regions[i].sub > 0 && regions[i].sub < subThreshold {
                regions[i].sub = 0
            }
        }

        return regions
    }

    // MARK: - Cross-Chain Coverage Helpers

    /// Maximum gap length at a given query length.
    /// Matches bwa-mem2's `cal_max_gap` (bwamem.cpp:66-76).
    private static func calMaxGap(queryLen: Int, scoring: ScoringParameters) -> Int {
        let a = Int(scoring.matchScore)
        let oDel = Int(scoring.gapOpenPenaltyDeletion)
        let eDel = Int(scoring.gapExtendPenaltyDeletion)
        let oIns = Int(scoring.gapOpenPenalty)
        let eIns = Int(scoring.gapExtendPenalty)
        let lDel = eDel > 0 ? (queryLen * a - oDel) / eDel + 1 : 1
        let lIns = eIns > 0 ? (queryLen * a - oIns) / eIns + 1 : 1
        var l = max(lDel, lIns)
        l = max(l, 1)
        return min(l, Int(scoring.bandWidth) << 1)
    }

    /// Check if a seed is covered by an existing alignment region using
    /// query+reference containment and diagonal proximity.
    /// Matches bwa-mem2's check in mem_chain2aln_across_reads_V2 (lines 2930-2955).
    private static func isSeedCoveredByRegion(
        seed: MemSeed,
        seedQEnd: Int32,
        seedREnd: Int64,
        region: MemAlnReg,
        readLength: Int32,
        scoring: ScoringParameters
    ) -> Bool {
        // 1. Seed must be fully contained on both query AND reference
        guard seed.qbeg >= region.qb && seedQEnd <= region.qe
           && seed.rbeg >= region.rb && seedREnd <= region.re else { return false }

        // 2. Seed length vs region's original seed length (bwa-mem2 line 2944)
        if seed.len - region.seedLen0 > Int32(Float(readLength) * 0.1) { return false }

        // 3. Diagonal proximity check from the front
        let qd = Int(seed.qbeg - region.qb)
        let rd = Int(seed.rbeg - region.rb)
        let minDist = min(qd, rd)
        var maxGap = calMaxGap(queryLen: minDist, scoring: scoring)
        var w = min(maxGap, Int(region.w))
        if abs(qd - rd) < w { return true }

        // 4. Diagonal proximity check from the back
        let qd2 = Int(region.qe - seedQEnd)
        let rd2 = Int(region.re - seedREnd)
        let minDist2 = min(qd2, rd2)
        maxGap = calMaxGap(queryLen: minDist2, scoring: scoring)
        w = min(maxGap, Int(region.w))
        if abs(qd2 - rd2) < w { return true }

        return false
    }

    /// Check if any other seed in the chain overlaps with the given seed on query
    /// coordinates but maps to a different reference diagonal, suggesting a genuinely
    /// different alignment. Matches bwa-mem2 lines 2965-2978.
    private static func hasOverlappingSeedEvidence(
        seedIdx: Int,
        seeds: SeedSoA,
        currentQbeg: Int32,
        currentLen: Int32,
        currentRbeg: Int64
    ) -> Bool {
        let q = seeds.qbegs, r = seeds.rbegs, l = seeds.lens
        for otherIdx in 0..<seeds.count {
            if otherIdx == seedIdx { continue }
            // Only check seeds of similar length
            if l[otherIdx] < Int32(Float(currentLen) * 0.95) { continue }
            // Check overlap from current seed's perspective
            if currentQbeg <= q[otherIdx]
                && currentQbeg + currentLen - q[otherIdx] >= currentLen >> 2
                && q[otherIdx] - currentQbeg != Int32(truncatingIfNeeded: r[otherIdx] - currentRbeg) {
                return true
            }
            // Check overlap from other seed's perspective
            if q[otherIdx] <= currentQbeg
                && q[otherIdx] + l[otherIdx] - currentQbeg >= currentLen >> 2
                && currentQbeg - q[otherIdx] != Int32(truncatingIfNeeded: currentRbeg - r[otherIdx]) {
                return true
            }
        }
        return false
    }

    // MARK: - GPU Batch Extension Support

    /// Collect extension plans for all seeds across all chains for a single read.
    /// Pre-fetches and reverses sequences for left/right extensions so they can be
    /// dispatched to the GPU in bulk. Does NOT perform covered-region filtering —
    /// all seeds are planned; filtering is applied post-hoc in `assembleRegions()`.
    public static func collectExtensionPlans(
        chains: [MemChain],
        read: ReadSequence,
        getReference: (Int64, Int) -> [UInt8],
        scoring: ScoringParameters
    ) -> [SeedExtensionPlan] {
        var plans: [SeedExtensionPlan] = []
        let readLen = Int32(read.length)
        let bandWidth = scoring.bandWidth

        for (chainIdx, chain) in chains.enumerated() {
            let seeds = chain.seeds
            for (seedIdx, seed) in seeds.enumerated() {
                let seedQEnd = seed.qbeg + seed.len
                let seedREnd = seed.rbeg + Int64(seed.len)

                // Left extension sequences (reversed)
                var leftQuery: [UInt8] = []
                var leftTarget: [UInt8] = []
                if seed.qbeg > 0 {
                    leftQuery = Array(read.bases[0..<Int(seed.qbeg)].reversed())
                    let targetLeftLen = min(Int(seed.qbeg) + Int(bandWidth), Int(seed.rbeg))
                    if targetLeftLen > 0 {
                        let targetStart = seed.rbeg - Int64(targetLeftLen)
                        leftTarget = Array(getReference(targetStart, targetLeftLen).reversed())
                    }
                }

                // Right extension sequences (forward)
                var rightQuery: [UInt8] = []
                var rightTarget: [UInt8] = []
                if seedQEnd < readLen {
                    rightQuery = Array(read.bases[Int(seedQEnd)..<read.length])
                    let targetRightLen = min(Int(readLen - seedQEnd) + Int(bandWidth), 10000)
                    if targetRightLen > 0 {
                        rightTarget = getReference(seedREnd, targetRightLen)
                    }
                }

                plans.append(SeedExtensionPlan(
                    chainIndex: chainIdx,
                    seedIndexInChain: seedIdx,
                    seed: seed,
                    seedQEnd: seedQEnd,
                    seedREnd: seedREnd,
                    readLen: readLen,
                    chainRid: chain.rid,
                    chainIsAlt: chain.isAlt,
                    chainFracRep: chain.fracRep,
                    chainSeedCount: seeds.count,
                    chainSeeds: seeds,
                    leftQuery: leftQuery,
                    leftTarget: leftTarget,
                    rightQuery: rightQuery,
                    rightTarget: rightTarget
                ))
            }
        }
        return plans
    }

    /// Assemble alignment regions from GPU extension results.
    /// Applies clip-vs-extend decisions, covered-region filtering, trueScore,
    /// and seedCov computation — exactly matching the CPU `extend()` logic.
    ///
    /// Plans must be ordered by (chainIndex, seedIndexInChain) — the order
    /// returned by `collectExtensionPlans()`.
    public static func assembleRegions(
        plans: [SeedExtensionPlan],
        leftResults: [SWResult],
        rightResults: [SWResult],
        scoring: ScoringParameters
    ) -> [MemAlnReg] {
        guard !plans.isEmpty else { return [] }

        let bandWidth = scoring.bandWidth
        // tmp = max single-event penalty (for subN threshold)
        let tmp = max(
            scoring.matchScore + scoring.mismatchPenalty,
            max(scoring.gapOpenPenalty + scoring.gapExtendPenalty,
                scoring.gapOpenPenaltyDeletion + scoring.gapExtendPenaltyDeletion)
        )

        var regions: [MemAlnReg] = []

        // Cache SeedSoA per unique chain index to avoid repeated conversion
        var chainSoACache: [Int: SeedSoA] = [:]
        defer { for (_, soa) in chainSoACache { soa.deallocate() } }

        for (planIdx, plan) in plans.enumerated() {
            let seed = plan.seed
            let readLen = plan.readLen

            // Get or create SeedSoA for this chain
            let chainSoA: SeedSoA
            if let cached = chainSoACache[plan.chainIndex] {
                chainSoA = cached
            } else {
                let soa = SeedSoA(from: plan.chainSeeds)
                chainSoACache[plan.chainIndex] = soa
                chainSoA = soa
            }

            // Cross-chain coverage check (matches CPU path)
            var coveringRegIdx = -1
            for (idx, reg) in regions.enumerated() {
                if isSeedCoveredByRegion(
                    seed: seed, seedQEnd: plan.seedQEnd, seedREnd: plan.seedREnd,
                    region: reg, readLength: readLen, scoring: scoring
                ) {
                    coveringRegIdx = idx
                    break
                }
            }

            if coveringRegIdx >= 0 {
                // Check for overlapping-seed evidence of different alignment
                if hasOverlappingSeedEvidence(
                    seedIdx: plan.seedIndexInChain, seeds: chainSoA,
                    currentQbeg: seed.qbeg, currentLen: seed.len, currentRbeg: seed.rbeg
                ) {
                    // Evidence found — keep this extension
                } else {
                    let seedAlnScore = seed.len * scoring.matchScore
                    if seedAlnScore > regions[coveringRegIdx].sub {
                        regions[coveringRegIdx].sub = seedAlnScore
                    }
                    if seedAlnScore >= regions[coveringRegIdx].score - tmp {
                        regions[coveringRegIdx].subN += 1
                    }
                    continue
                }
            }

            var reg = MemAlnReg()
            reg.w = bandWidth
            reg.rid = plan.chainRid
            reg.isAlt = plan.chainIsAlt
            reg.fracRep = plan.chainFracRep
            reg.seedCov = seed.len
            reg.seedLen0 = seed.len

            var accumulatedH0 = seed.score

            // --- Left extension result ---
            var leftScore: Int32 = 0
            var leftQLen: Int32 = 0
            var leftTLen: Int32 = 0

            if !plan.leftQuery.isEmpty && !plan.leftTarget.isEmpty {
                let result = leftResults[planIdx]
                accumulatedH0 = result.score

                if result.globalScore <= 0
                    || result.globalScore <= result.score - scoring.penClip5 {
                    leftQLen = result.queryEnd
                    leftTLen = result.targetEnd
                    leftScore = result.score
                } else {
                    leftQLen = seed.qbeg
                    leftTLen = result.globalTargetEnd
                    leftScore = result.globalScore
                }
            }

            // --- Right extension result ---
            var rightScore: Int32 = 0
            var rightQLen: Int32 = 0
            var rightTLen: Int32 = 0

            if !plan.rightQuery.isEmpty && !plan.rightTarget.isEmpty {
                let result = rightResults[planIdx]

                if result.globalScore <= 0
                    || result.globalScore <= result.score - scoring.penClip3 {
                    rightQLen = result.queryEnd
                    rightTLen = result.targetEnd
                    rightScore = result.score
                } else {
                    rightQLen = plan.readLen - plan.seedQEnd
                    rightTLen = result.globalTargetEnd
                    rightScore = result.globalScore
                }
            }

            // Compute final alignment coordinates
            reg.qb = seed.qbeg - leftQLen
            reg.qe = plan.seedQEnd + rightQLen
            reg.rb = seed.rbeg - Int64(leftTLen)
            reg.re = plan.seedREnd + Int64(rightTLen)

            var trueScore = Int32(seed.len) * scoring.matchScore
            if leftScore > 0 {
                trueScore = leftScore
            }
            if rightScore > 0 {
                trueScore += rightScore - accumulatedH0
            }
            reg.score = trueScore
            reg.trueScore = trueScore

            reg.seedCov = computeSeedCov(seeds: chainSoA, qb: reg.qb, qe: reg.qe, rb: reg.rb, re: reg.re)

            regions.append(reg)
        }

        // Threshold: sub must exceed minSeedLength * matchScore
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
