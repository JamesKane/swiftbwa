import BWACore
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif canImport(Musl)
import Musl
#endif

/// Result of pair resolution for one read pair.
public struct PairDecision: Sendable {
    /// Index into regions1
    public var idx1: Int
    /// Index into regions2
    public var idx2: Int
    /// Combined pair score (best pair)
    public var score: Int
    /// Second-best pair score (0 if only one valid pairing)
    public var secondBestScore: Int
    /// Count of pairs within score threshold of second-best
    public var nSub: Int
    /// Whether this pair is proper (concordant orientation + reasonable insert size)
    public var isProperPair: Bool
    /// Insert size between the pair
    public var insertSize: Int64
    /// Pair orientation
    public var orientation: PairOrientation

    public init(idx1: Int, idx2: Int, score: Int, secondBestScore: Int = 0,
                nSub: Int = 0, isProperPair: Bool,
                insertSize: Int64, orientation: PairOrientation) {
        self.idx1 = idx1
        self.idx2 = idx2
        self.score = score
        self.secondBestScore = secondBestScore
        self.nSub = nSub
        self.isProperPair = isProperPair
        self.insertSize = insertSize
        self.orientation = orientation
    }
}

/// Resolves paired-end alignments using z-score-based pair scoring.
///
/// Reimplements `mem_pair()` from bwa-mem2's `bwamem_pair.cpp`.
public struct PairedEndResolver: Sendable {

    /// Z-score threshold for proper pairing.
    private static let zThreshold: Double = 4.0

    /// Resolve the best pair from two sets of alignment regions.
    ///
    /// Algorithm (from `mem_pair()` in bwamem_pair.cpp:277-370):
    /// 1. Consider top `maxPairings` primary/sub-primary regions per read
    /// 2. Score each (r1, r2) combination using z-score-based penalty
    /// 3. Pick pair with highest score
    ///
    /// - Parameters:
    ///   - regions1: Alignment regions for read 1
    ///   - regions2: Alignment regions for read 2
    ///   - dist: Insert size distribution (from InsertSizeEstimator)
    ///   - genomeLength: Forward genome length
    ///   - scoring: Scoring parameters
    /// - Returns: Best pair decision, or nil if no valid pairing exists
    public static func resolve(
        regions1: [MemAlnReg],
        regions2: [MemAlnReg],
        dist: InsertSizeDistribution,
        genomeLength: Int64,
        scoring: ScoringParameters
    ) -> PairDecision? {
        guard !regions1.isEmpty && !regions2.isEmpty else { return nil }

        // Select pairing candidates (primaries + secondaries within score threshold)
        let candidates1 = topCandidates(from: regions1, scoring: scoring)
        let candidates2 = topCandidates(from: regions2, scoring: scoring)

        // Collect all valid pair scores (bwa-mem2 mem_pair u[] array)
        struct PairHit {
            var score: Int
            var idx1: Int
            var idx2: Int
            var isProper: Bool
            var insertSize: Int64
            var orientation: PairOrientation
        }
        var hits: [PairHit] = []

        for (i, origIdx1) in candidates1 {
            let r1 = regions1[i]

            for (j, origIdx2) in candidates2 {
                let r2 = regions2[j]

                guard let result = InsertSizeEstimator.inferOrientation(
                    r1: r1, r2: r2, genomeLength: genomeLength
                ) else { continue }

                let oriStats = dist.stats[result.orientation.rawValue]

                // Only consider pairs within the proper-pair distance window
                guard !oriStats.failed && oriStats.stddev > 0
                    && result.insertSize >= oriStats.properLow
                    && result.insertSize <= oriStats.properHigh else { continue }

                let baseScore = Int(r1.score) + Int(r2.score)
                let z = abs(Double(result.insertSize) - oriStats.mean) / oriStats.stddev
                let penalty = computeZPenalty(z: z, matchScore: scoring.matchScore)
                let pairScore = baseScore - penalty
                let isProper = (result.orientation == dist.primaryOrientation)

                hits.append(PairHit(
                    score: pairScore, idx1: origIdx1, idx2: origIdx2,
                    isProper: isProper, insertSize: result.insertSize,
                    orientation: result.orientation
                ))
            }
        }

        guard !hits.isEmpty else { return nil }

        // Sort descending by score to find best and second-best
        hits.sort { $0.score > $1.score }

        let best = hits[0]

        // Second-best score and n_sub (bwa-mem2 mem_pair lines 337-342)
        let subo = hits.count > 1 ? hits[1].score : 0
        let tmp = Int(scoring.matchScore) + Int(scoring.mismatchPenalty)
        var nSub = 0
        for i in 1..<hits.count {
            if subo - hits[i].score <= tmp {
                nSub += 1
            }
        }

        return PairDecision(
            idx1: best.idx1,
            idx2: best.idx2,
            score: best.score,
            secondBestScore: subo,
            nSub: nSub,
            isProperPair: best.isProper,
            insertSize: best.insertSize,
            orientation: best.orientation
        )
    }

    /// Compute TLEN values for a read pair per SAM spec.
    ///
    /// - Parameters:
    ///   - pos1: 0-based reference position of read 1
    ///   - isReverse1: Whether read 1 is on reverse strand
    ///   - refLen1: Reference bases consumed by read 1's CIGAR
    ///   - pos2: 0-based reference position of read 2
    ///   - isReverse2: Whether read 2 is on reverse strand
    ///   - refLen2: Reference bases consumed by read 2's CIGAR
    /// - Returns: TLEN values for read 1 and read 2
    public static func computeTLEN(
        pos1: Int64, isReverse1: Bool, refLen1: Int64,
        pos2: Int64, isReverse2: Bool, refLen2: Int64
    ) -> (tlen1: Int64, tlen2: Int64) {
        let end1 = pos1 + refLen1
        let end2 = pos2 + refLen2

        let leftmost = min(pos1, pos2)
        let rightmost = max(end1, end2)
        let span = rightmost - leftmost

        if pos1 < pos2 || (pos1 == pos2 && !isReverse1) {
            // Read 1 is leftmost
            return (span, -span)
        } else {
            // Read 2 is leftmost
            return (-span, span)
        }
    }

    /// Adjust per-end MAPQ using paired-end evidence.
    ///
    /// Reimplements bwa-mem2's PE MAPQ logic from `mem_sam_pe()` (bwamem_pair.cpp:441-467).
    ///
    /// - Parameters:
    ///   - seMAPQ1: Single-end MAPQ for read 1
    ///   - seMAPQ2: Single-end MAPQ for read 2
    ///   - region1: Best alignment region for read 1 (after pair promotion)
    ///   - region2: Best alignment region for read 2 (after pair promotion)
    ///   - decision: Pair resolution result (best pair score, subo, n_sub)
    ///   - scoring: Scoring parameters
    /// - Returns: Adjusted (mapq1, mapq2)
    public static func adjustMAPQ(
        seMAPQ1: UInt8,
        seMAPQ2: UInt8,
        region1: MemAlnReg,
        region2: MemAlnReg,
        decision: PairDecision,
        scoring: ScoringParameters
    ) -> (UInt8, UInt8) {
        let o = decision.score
        var subo = decision.secondBestScore

        // bwa-mem2 line 441-443: compare against unpaired SE score
        let scoreUn = Int(region1.score) + Int(region2.score) - Int(scoring.unpairedPenalty)
        if subo < scoreUn { subo = scoreUn }

        // bwa-mem2 line 444: q_pe = raw_mapq(o - subo, a)
        let a = Double(scoring.matchScore)
        var qPe = Int(6.02 * Double(o - subo) / a + 0.499)

        // bwa-mem2 line 445-446: n_sub penalty
        if decision.nSub > 0 {
            qPe -= Int(4.343 * log(Double(decision.nSub + 1)) + 0.499)
        }

        // Clamp q_pe
        if qPe < 0 { qPe = 0 }
        if qPe > 60 { qPe = 60 }

        // frac_rep not implemented (would be: qPe *= (1 - 0.5*(frac_rep1 + frac_rep2)))

        var q1 = Int(seMAPQ1)
        var q2 = Int(seMAPQ2)

        if o > scoreUn {
            // bwa-mem2 lines 461-462: paired alignment preferred — boost SE MAPQ with pair evidence
            // q_se = max(q_se, min(q_pe, q_se + 40))
            q1 = q1 > qPe ? q1 : (qPe < q1 + 40 ? qPe : q1 + 40)
            q2 = q2 > qPe ? q2 : (qPe < q2 + 40 ? qPe : q2 + 40)

            // bwa-mem2 lines 466-467: cap at tandem repeat score (csub)
            if region1.csub > 0 {
                let csubCap = Int(6.02 * Double(Int(region1.score) - Int(region1.csub)) / a + 0.499)
                if q1 > csubCap { q1 = csubCap }
            }
            if region2.csub > 0 {
                let csubCap = Int(6.02 * Double(Int(region2.score) - Int(region2.csub)) / a + 0.499)
                if q2 > csubCap { q2 = csubCap }
            }
        }
        // else: unpaired preferred — keep SE MAPQ as-is (bwa-mem2 lines 469-472)

        // Final clamp
        q1 = max(0, min(60, q1))
        q2 = max(0, min(60, q2))

        return (UInt8(q1), UInt8(q2))
    }

    // MARK: - Private Helpers

    /// Select pairing candidates: all primary regions, plus secondaries within
    /// score threshold. For multi-mappers (MAPQ=0), all positions have the same
    /// score but most are secondary; including them allows proper-pair discovery.
    /// Returns (index_in_input_array, original_index) pairs.
    private static func topCandidates(
        from regions: [MemAlnReg],
        scoring: ScoringParameters
    ) -> [(Int, Int)] {
        guard let bestScore = regions.first(where: { $0.secondary < 0 })?.score else {
            // No primary — use all regions
            return regions.indices.map { ($0, $0) }
        }
        let threshold = bestScore - scoring.unpairedPenalty
        var result: [(Int, Int)] = []
        for (i, region) in regions.enumerated() {
            guard region.score >= threshold else { continue }
            result.append((i, i))
        }
        return result
    }

    /// Compute z-score-based penalty using bwa-mem2's exact formula:
    /// `score1 + score2 + 0.721 * ln(2 * erfc(|z| / sqrt(2))) * matchScore + 0.499`
    /// where 0.721 ≈ 1/ln(4). The log term is negative, reducing the pair score.
    private static func computeZPenalty(z: Double, matchScore: Int32) -> Int {
        let erfcVal = erfc(z / sqrt(2.0))
        guard erfcVal > 0 else {
            // z is so large that erfc underflows to 0; treat as maximum penalty
            return Int(matchScore) * 10
        }
        let logVal = log(2.0 * erfcVal)
        let penalty = 0.721 * logVal * Double(matchScore) + 0.499
        // Penalty should be non-positive (it reduces score); floor to int
        return max(0, Int(floor(-penalty)))
    }
}
