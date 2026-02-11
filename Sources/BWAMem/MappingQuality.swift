import BWACore
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif canImport(Musl)
import Musl
#endif

/// Computes mapping quality (MAPQ) from alignment scores.
/// Reimplements BWA-MEM's MAPQ calculation logic.
public struct MappingQuality: Sendable {

    /// Compute MAPQ for a primary alignment.
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

        let maxScore = readLength * scoring.matchScore

        // If unmapped or zero score
        if region.score <= 0 { return 0 }

        // Perfect unique hit
        if allRegions.count == 1 && region.score == maxScore {
            return 60
        }

        // Calculate identity-based component
        let identity = Float(region.score) / Float(maxScore)
        var mapq: Float = 0

        if allRegions.count <= 1 || region.secondary < 0 {
            // Primary alignment
            let subScore = region.sub > 0 ? region.sub : 0
            let scoreDiff = region.score - subScore

            if scoreDiff == 0 {
                mapq = 0
            } else {
                // BWA formula: mapq ~ 250 * (1 - sub/score) * min(1, score/maxScore) * log2(score - sub)
                let ratio = 1.0 - Float(subScore) / Float(region.score)
                let logComp = log2f(Float(scoreDiff) + 1.0)
                mapq = 250.0 * ratio * identity * logComp / log2f(Float(maxScore) + 1.0)
            }

            // Cap at 60
            mapq = min(mapq, 60)

            // Penalize if sub_n > 0
            if region.subN > 0 {
                mapq -= 4.343 * log(1.0 + Float(region.subN))
            }
        } else {
            // Secondary alignment - always low MAPQ
            mapq = 0
        }

        return UInt8(clamping: max(0, Int(mapq)))
    }
}
