import BWACore
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif canImport(Musl)
import Musl
#endif

/// bwa-mem2 defaults: mapQ_coef_len = 50, mapQ_coef_fac = (int)log(50) = 3
/// Note: bwa-mem2 declares mapQ_coef_fac as int, truncating log(50)=3.912 to 3.
private let mapQCoefLen: Double = 50.0
private let mapQCoefFac: Double = 3.0

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
        // bwa-mem2 line 1474: tandem repeat detection — csub from mate rescue score2
        if region.csub > sub { sub = region.csub }

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
            // bwa-mem2 line 1480-1484 (mapQ_coef_len=50 > 0, default branch):
            // tmp = l < mapQ_coef_len ? 1. : mapQ_coef_fac / log(l);
            // tmp *= identity * identity;
            // mapq = (int)(6.02 * (score - sub) / a * tmp * tmp + .499);
            let dl = Double(l)
            var tmp = dl < mapQCoefLen ? 1.0 : mapQCoefFac / log(dl)
            tmp *= identity * identity
            mapq = Int32(6.02 * Double(region.score - sub) / a * tmp * tmp + 0.499)
        }

        // bwa-mem2 line 1489: sub_n penalty
        if region.subN > 0 {
            mapq -= Int32(4.343 * log(Double(region.subN + 1)) + 0.499)
        }

        // bwa-mem2 lines 1490-1491: clamp to [0, 60]
        if mapq > 60 { mapq = 60 }
        if mapq < 0 { mapq = 0 }

        // bwa-mem2 line 1492: frac_rep scaling
        mapq = Int32(Double(mapq) * (1.0 - Double(region.fracRep)) + 0.499)

        return UInt8(clamping: mapq)
    }
}
