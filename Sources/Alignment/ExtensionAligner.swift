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
        let mat = scoring.scoringMatrix()
        let queryBuf = UnsafeMutablePointer<UInt8>.allocate(capacity: read.length)
        let targetBuf = UnsafeMutablePointer<UInt8>.allocate(capacity: 10200)
        defer { queryBuf.deallocate(); targetBuf.deallocate() }
        return extend(
            chain: chain, read: read,
            getReference: { pos, length, buf in
                let bases = getReference(pos, length)
                bases.withUnsafeBufferPointer { src in
                    buf.initialize(from: src.baseAddress!, count: src.count)
                }
                return bases.count
            },
            scoring: scoring, scoringMatrix: mat,
            queryBuf: queryBuf, targetBuf: targetBuf
        )
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
        let lQuery = Int(readLen)

        // Compute chain-level reference range (bwa-mem2 bwamem.cpp:2146-2157)
        var rmax0: Int64 = Int64.max
        var rmax1: Int64 = 0
        for i in 0..<soa.count {
            let t = (qbeg: soa.qbegs[i], len: soa.lens[i], rbeg: soa.rbegs[i])
            let b = t.rbeg - Int64(Int(t.qbeg) + calMaxGap(queryLen: Int(t.qbeg), scoring: scoring))
            let remaining = lQuery - Int(t.qbeg) - Int(t.len)
            let e = t.rbeg + Int64(t.len) + Int64(remaining + calMaxGap(queryLen: remaining, scoring: scoring))
            if b < rmax0 { rmax0 = b }
            if e > rmax1 { rmax1 = e }
        }
        if rmax0 < 0 { rmax0 = 0 }

        let tmp = max(
            scoring.matchScore + scoring.mismatchPenalty,
            max(scoring.gapOpenPenalty + scoring.gapExtendPenalty,
                scoring.gapOpenPenaltyDeletion + scoring.gapExtendPenaltyDeletion)
        )

        // Sort seeds by score descending (bwa-mem2 bwamem.cpp:2189-2207)
        // Highest-scoring (longest) seeds get extended first, matching bwa-mem2's
        // `srt[i] = score<<32 | i` with reverse iteration.
        var seedOrder = Array(0..<soa.count)
        if soa.count > 1 {
            seedOrder.sort { soa.scores[$0] > soa.scores[$1] }
        }

        for sortedIdx in seedOrder {
            let seedQbeg = soa.qbegs[sortedIdx]
            let seedLen = soa.lens[sortedIdx]
            let seedRbeg = soa.rbegs[sortedIdx]
            let seedScore = soa.scores[sortedIdx]

            let seedQEnd = seedQbeg + seedLen
            let seedREnd = seedRbeg + Int64(seedLen)

            // Covered-region check: query containment
            // Seed ordering by score (above) ensures longest (best-scoring) seeds
            // get extended first, matching bwa-mem2's seed processing order.
            var coveringRegIdx = -1
            for (idx, reg) in regions.enumerated() {
                if seedQbeg >= reg.qb && seedQEnd <= reg.qe {
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

            // --- Left extension (with band retry, MAX_BAND_TRY=2) ---
            var leftScore: Int32 = 0
            var leftQLen: Int32 = 0
            var leftTLen: Int32 = 0

            if seedQbeg > 0 {
                let queryLeftLen = Int(seedQbeg)
                // Write reversed query into queryBuf
                for i in 0..<queryLeftLen {
                    queryBuf[i] = read.bases[queryLeftLen - 1 - i]
                }

                // Target length from chain rmax (bwa-mem2 bwamem.cpp:2281)
                let targetLeftLen = min(Int(seedRbeg - rmax0), 10000)
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

                    // Band retry loop (bwa-mem2 bwamem.cpp:2473-2527)
                    var w = bandWidth
                    var prevScore = seedScore
                    var leftResult = SWResult()
                    for tryIdx: Int32 in 0..<2 {
                        leftResult = BandedSWScalar.align(
                            query: qBuf, target: tBuf,
                            scoring: scoring, w: w, h0: seedScore,
                            endBonus: scoring.penClip5
                        )
                        if leftResult.score == prevScore
                            || leftResult.maxOff < (w >> 1) + (w >> 2)
                            || tryIdx == 1 {
                            break
                        }
                        prevScore = leftResult.score
                        w = w << 1
                    }

                    accumulatedH0 = leftResult.score
                    reg.w = max(reg.w, w)

                    if leftResult.globalScore <= 0
                        || leftResult.globalScore <= leftResult.score - scoring.penClip5 {
                        leftQLen = leftResult.queryEnd
                        leftTLen = leftResult.targetEnd
                        leftScore = leftResult.score
                    } else {
                        leftQLen = seedQbeg
                        leftTLen = leftResult.globalTargetEnd
                        leftScore = leftResult.globalScore
                    }
                }
            }

            // --- Right extension (with band retry, MAX_BAND_TRY=2) ---
            var rightScore: Int32 = 0
            var rightQLen: Int32 = 0
            var rightTLen: Int32 = 0

            if seedQEnd < readLen {
                let queryRightLen = read.length - Int(seedQEnd)
                // Write forward query into queryBuf
                for i in 0..<queryRightLen {
                    queryBuf[i] = read.bases[Int(seedQEnd) + i]
                }

                // Target length from chain rmax (bwa-mem2 bwamem.cpp:2357)
                let targetRightLen = min(Int(rmax1 - seedREnd), 10000)
                if targetRightLen > 0 {
                    let actualTargetLen = getReference(seedREnd, targetRightLen, targetBuf)

                    let qBuf = UnsafeBufferPointer(start: queryBuf, count: queryRightLen)
                    let tBuf = UnsafeBufferPointer(start: targetBuf, count: actualTargetLen)

                    // Band retry loop (bwa-mem2 bwamem.cpp:2605-2660)
                    var w = bandWidth
                    var prevScore = accumulatedH0
                    var rightResult = SWResult()
                    for tryIdx: Int32 in 0..<2 {
                        rightResult = BandedSWScalar.align(
                            query: qBuf, target: tBuf,
                            scoring: scoring, w: w, h0: accumulatedH0,
                            endBonus: scoring.penClip3
                        )
                        if rightResult.score == prevScore
                            || rightResult.maxOff < (w >> 1) + (w >> 2)
                            || tryIdx == 1 {
                            break
                        }
                        prevScore = rightResult.score
                        w = w << 1
                    }

                    reg.w = max(reg.w, w)

                    if rightResult.globalScore <= 0
                        || rightResult.globalScore <= rightResult.score - scoring.penClip3 {
                        rightQLen = rightResult.queryEnd
                        rightTLen = rightResult.targetEnd
                        rightScore = rightResult.score
                    } else {
                        rightQLen = readLen - seedQEnd
                        rightTLen = rightResult.globalTargetEnd
                        rightScore = rightResult.globalScore
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

}
