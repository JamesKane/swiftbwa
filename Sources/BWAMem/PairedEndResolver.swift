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
    /// Combined pair score
    public var score: Int
    /// Whether this pair is proper (concordant orientation + reasonable insert size)
    public var isProperPair: Bool
    /// Insert size between the pair
    public var insertSize: Int64
    /// Pair orientation
    public var orientation: PairOrientation

    public init(idx1: Int, idx2: Int, score: Int, isProperPair: Bool,
                insertSize: Int64, orientation: PairOrientation) {
        self.idx1 = idx1
        self.idx2 = idx2
        self.score = score
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

        var bestDecision: PairDecision?
        var bestScore = Int.min

        for (i, origIdx1) in candidates1 {
            let r1 = regions1[i]

            for (j, origIdx2) in candidates2 {
                let r2 = regions2[j]

                guard let result = InsertSizeEstimator.inferOrientation(
                    r1: r1, r2: r2, genomeLength: genomeLength
                ) else { continue }

                let oriStats = dist.stats[result.orientation.rawValue]

                // Only consider pairs within the proper-pair distance window
                // (bwa-mem2 mem_pair lines 318-319: dist < low → skip, dist > high → break)
                guard !oriStats.failed && oriStats.stddev > 0
                    && result.insertSize >= oriStats.properLow
                    && result.insertSize <= oriStats.properHigh else { continue }

                let baseScore = Int(r1.score) + Int(r2.score)
                let z = abs(Double(result.insertSize) - oriStats.mean) / oriStats.stddev
                let penalty = computeZPenalty(z: z, matchScore: scoring.matchScore)
                let pairScore = baseScore - penalty
                let isProper = (result.orientation == dist.primaryOrientation)

                if pairScore > bestScore {
                    bestScore = pairScore
                    bestDecision = PairDecision(
                        idx1: origIdx1,
                        idx2: origIdx2,
                        score: pairScore,
                        isProperPair: isProper,
                        insertSize: result.insertSize,
                        orientation: result.orientation
                    )
                }
            }
        }

        return bestDecision
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

    /// Adjust MAPQ for paired-end reads.
    ///
    /// From bwa-mem2 `mem_sam_pe()`: for proper pairs, MAPQ is boosted if the
    /// pair score gap implies high confidence.
    ///
    /// - Parameters:
    ///   - mapq: Original single-end MAPQ
    ///   - pairScore: Score of the best pair
    ///   - secondBestPairScore: Score of the second-best pair (nil if only one pairing)
    ///   - isProperPair: Whether the best pair is proper
    /// - Returns: Adjusted MAPQ
    public static func adjustMAPQ(
        mapq: UInt8,
        pairScore: Int,
        secondBestPairScore: Int?,
        isProperPair: Bool
    ) -> UInt8 {
        guard isProperPair else { return mapq }

        if let secondBest = secondBestPairScore {
            let gap = pairScore - secondBest
            if gap > 0 {
                let pairMAPQ = min(60, Int(4.343 * log(Double(gap)) + 0.5))
                return UInt8(min(60, max(Int(mapq), pairMAPQ)))
            }
        } else {
            // No second-best pair: boost to at least 30 for proper pairs
            return UInt8(min(60, max(Int(mapq), 30)))
        }
        return mapq
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
