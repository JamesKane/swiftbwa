import BWACore
import FMIndex
import Alignment

/// Mate rescue: Smith-Waterman alignment of an unmapped mate near the expected
/// position based on the mapped read and insert size distribution.
///
/// Equivalent to bwa-mem2's `mem_matesw()` from `bwamem_pair.cpp:150-283`.
public struct MateRescue: Sendable {

    /// Attempt to rescue alignments for `mateRead` using `templateRegions` as anchors.
    ///
    /// For each template alignment and each of 4 bwa directions, computes the expected
    /// rescue region, extracts reference, runs local SW, and returns any rescued regions.
    ///
    /// - Parameters:
    ///   - templateRegions: Filtered regions from the mapped read (score >= best - penUnpaired)
    ///   - mateRead: The read to rescue
    ///   - mateRegions: Existing alignment regions for the mate (to check for redundancy)
    ///   - dist: Insert size distribution
    ///   - genomeLength: Forward genome length
    ///   - packedRef: Packed reference for sequence extraction
    ///   - metadata: Reference metadata for chromosome boundaries
    ///   - scoring: Scoring parameters
    /// - Returns: Array of rescued alignment regions
    public static func rescue(
        templateRegions: [MemAlnReg],
        mateRead: ReadSequence,
        mateRegions: [MemAlnReg],
        dist: InsertSizeDistribution,
        genomeLength: Int64,
        packedRef: PackedReference,
        metadata: ReferenceMetadata,
        scoring: ScoringParameters,
        scoringMatrix: [Int8]? = nil
    ) -> [MemAlnReg] {
        guard !templateRegions.isEmpty else { return [] }

        let primaryStats = dist.stats[dist.primaryOrientation.rawValue]
        guard !primaryStats.failed else { return [] }

        let mateLen = Int64(mateRead.length)
        let minScore = scoring.minSeedLength * scoring.matchScore
        var rescued: [MemAlnReg] = []

        for a in templateRegions {
            // For each of 4 bwa directions (FF=0, FR=1, RF=2, RR=3)
            var allDirectionsSatisfied = true

            for r in 0..<4 {
                let oriStats = dist.stats[r]
                // Skip orientations with insufficient or unreliable data.
                // Non-primary orientations with few observations cause
                // excessive SW calls (O(N) per pair) with negligible benefit.
                guard !oriStats.failed else { continue }

                // Check if mate already has a hit for this direction
                if hasMateHit(mateRegions: mateRegions, template: a,
                              direction: r, stats: oriStats, genomeLength: genomeLength) {
                    continue
                }
                allDirectionsSatisfied = false

                // Compute rescue region
                // is_rev: mate is on opposite strand from template
                let isRev = ((r >> 1) != (r & 1))
                // is_larger: mate has larger BWT coordinate
                let isLarger = ((r >> 1) == 0)

                let pLow = oriStats.properLow
                let pHigh = oriStats.properHigh

                var rb: Int64
                var re: Int64

                if !isRev {
                    // Same strand rescue
                    if isLarger {
                        rb = a.rb + pLow
                        re = a.rb + pHigh + mateLen
                    } else {
                        rb = a.rb - pHigh - mateLen
                        re = a.rb - pLow
                    }
                } else {
                    // Opposite strand rescue
                    if isLarger {
                        rb = a.rb + pLow - mateLen
                        re = a.rb + pHigh
                    } else {
                        rb = a.rb - pHigh
                        re = a.rb - pLow + mateLen
                    }
                }

                // Clamp to valid BWT range
                rb = max(rb, 0)
                re = min(re, 2 * genomeLength)
                guard rb < re else { continue }

                // Extract reference sequence (strand-aware)
                guard let fetched = packedRef.fetchSequence(
                    bwtBegin: rb, bwtEnd: re,
                    genomeLength: genomeLength, metadata: metadata
                ) else { continue }

                // Must be same chromosome as template
                guard fetched.rid == a.rid else { continue }

                let refBases = fetched.bases
                let clampedRb = fetched.clampedBegin

                // Prepare mate sequence
                let mateSeq: [UInt8]
                if isRev {
                    mateSeq = mateRead.reverseComplement()
                } else {
                    mateSeq = mateRead.bases
                }

                // Run local SW
                guard let swResult = LocalSWAligner.align(
                    query: mateSeq, target: refBases, scoring: scoring,
                    scoringMatrix: scoringMatrix
                ) else { continue }

                // Filter by minimum score
                guard swResult.score >= minScore else { continue }

                // Construct MemAlnReg with coordinate conversion
                var reg = MemAlnReg()
                reg.score = swResult.score
                reg.trueScore = swResult.score
                reg.rid = a.rid
                reg.secondary = -1

                let tb = Int64(swResult.targetBegin)
                let te = Int64(swResult.targetEnd)
                let qb = swResult.queryBegin
                let qe = swResult.queryEnd

                if !isRev {
                    // Forward: direct coordinate mapping
                    reg.rb = clampedRb + tb
                    reg.re = clampedRb + te + 1
                    reg.qb = qb
                    reg.qe = qe + 1
                } else {
                    // Reverse: convert back to BWT reverse-strand coordinates
                    reg.rb = 2 * genomeLength - (clampedRb + te + 1)
                    reg.re = 2 * genomeLength - (clampedRb + tb)
                    reg.qb = Int32(mateLen) - (qe + 1)
                    reg.qe = Int32(mateLen) - qb
                }

                // Approximate seed coverage
                let refSpan = te - tb + 1
                let querySpan = Int64(qe - qb + 1)
                reg.seedCov = Int32(min(refSpan, querySpan) / 2)

                rescued.append(reg)
            }

            if allDirectionsSatisfied { break }
        }

        return rescued
    }

    // MARK: - Private Helpers

    /// Check if existing mate regions already have a hit consistent with the given
    /// direction and insert size bounds.
    private static func hasMateHit(
        mateRegions: [MemAlnReg],
        template: MemAlnReg,
        direction: Int,
        stats: OrientationStats,
        genomeLength: Int64
    ) -> Bool {
        for m in mateRegions {
            guard m.rid == template.rid else { continue }

            // Compute the bwa direction between template and this mate region
            let sameStrand = (template.rb >= genomeLength) == (m.rb >= genomeLength)
            let p2 = sameStrand ? m.rb : (2 * genomeLength - 1 - m.rb)
            let bwaDir = (sameStrand ? 0 : 1) ^ (p2 > template.rb ? 0 : 3)

            guard bwaDir == direction else { continue }

            let dist = abs(p2 - template.rb)
            if dist >= stats.properLow && dist <= stats.properHigh {
                return true
            }
        }
        return false
    }
}
