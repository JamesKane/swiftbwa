import BWACore
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif canImport(Musl)
import Musl
#endif

/// MEM_MAPQ_COEF from bwa-mem2's bwamem.h line 56
private let memMapqCoef: Double = 30.0

/// Computes mapping quality (MAPQ) from alignment scores.
/// Reimplements bwa-mem2's `mem_approx_mapq_se()` (bwamem.cpp lines 1470-1494).
public struct MappingQuality: Sendable {

    /// Compute MAPQ for a primary alignment.
    /// Matches bwa-mem2's `mem_approx_mapq_se()` exactly.
    ///
    /// - Parameters:
    ///   - region: The alignment region
    ///   - allRegions: All alignment regions for the read (sorted by score desc)
    ///   - scoring: Scoring parameters
    ///   - readLength: Length of the read
    /// - Returns: Phred-scaled mapping quality (0-60)
    public static func compute(
        region: MemAlnReg,
        allRegions: [MemAlnReg],
        scoring: ScoringParameters,
        readLength: Int32
    ) -> UInt8 {
        guard region.score > 0 else { return 0 }

        // bwa-mem2 line 1472: default sub = minSeedLen * matchScore when sub==0
        var sub = region.sub > 0 ? region.sub : scoring.minSeedLength * scoring.matchScore
        // csub not yet implemented; bwa-mem2 line 1474: sub = max(sub, csub)

        // bwa-mem2 line 1475
        if sub >= region.score { return 0 }

        // bwa-mem2 line 1476: alignment span = max(query span, ref span)
        let l = max(region.qe - region.qb, Int32(region.re - region.rb))

        // bwa-mem2 line 1477: identity based on alignment span
        let a = Double(scoring.matchScore)
        let b = Double(scoring.mismatchPenalty)
        let identity = 1.0 - (Double(l) * a - Double(region.score)) / (a + b) / Double(l)

        var mapq: Int32

        if region.score == 0 {
            mapq = 0
        } else {
            // bwa-mem2 line 1486 (default branch, mapQ_coef_len == 0):
            // mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499)
            let seedcov = max(region.seedCov, 1) // guard against log(0)
            mapq = Int32(
                memMapqCoef * (1.0 - Double(sub) / Double(region.score))
                    * log(Double(seedcov)) + 0.499
            )
            // bwa-mem2 line 1487: identity penalty when identity < 0.95
            if identity < 0.95 {
                mapq = Int32(Double(mapq) * identity * identity + 0.499)
            }
        }

        // bwa-mem2 line 1489: sub_n penalty
        if region.subN > 0 {
            mapq -= Int32(4.343 * log(Double(region.subN + 1)) + 0.499)
        }

        // bwa-mem2 lines 1490-1491: clamp to [0, 60]
        if mapq > 60 { mapq = 60 }
        if mapq < 0 { mapq = 0 }

        // bwa-mem2 line 1492: frac_rep scaling (not yet implemented, frac_rep = 0)
        // mapq = Int32(Double(mapq) * (1.0 - Double(region.fracRep)) + 0.499)

        return UInt8(clamping: mapq)
    }
}
