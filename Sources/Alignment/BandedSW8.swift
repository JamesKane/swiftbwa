import BWACore

/// 8-bit SIMD banded Smith-Waterman using SIMD16<UInt8> (128-bit NEON).
/// Processes 16 anti-diagonal cells simultaneously.
/// Falls back to 16-bit or scalar if scores overflow UInt8.
public struct BandedSW8: Sendable {

    /// Perform 8-bit SIMD banded Smith-Waterman.
    /// Returns nil if scores overflow 8-bit range (caller should use BandedSW16 or scalar).
    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        h0: Int32
    ) -> SWResult? {
        let qlen = query.count
        let tlen = target.count

        guard qlen > 0 && tlen > 0 else { return SWResult() }

        // Check if scores can fit in 8-bit
        let maxPossibleScore = h0 + Int32(qlen) * scoring.matchScore
        if maxPossibleScore > 250 {
            return nil  // Overflow risk, use 16-bit
        }

        let mat = scoring.scoringMatrix()
        let m = 5  // alphabet size
        let lanes = 16  // SIMD16<UInt8> width

        // Number of SIMD stripes needed to cover query
        let stripeCount = (qlen + lanes - 1) / lanes

        guard stripeCount > 0 else { return SWResult() }

        let zeroVec = SIMD16<UInt8>(repeating: 0)

        // Build striped query profile
        // For each target base k, for each stripe s, build SIMD16 of scores
        var profile = [[SIMD16<UInt8>]](repeating: [SIMD16<UInt8>](repeating: zeroVec, count: stripeCount), count: m)

        let bias = UInt8(scoring.mismatchPenalty)  // Add bias to keep values unsigned

        for k in 0..<m {
            for s in 0..<stripeCount {
                var vec = SIMD16<UInt8>(repeating: 0)
                for lane in 0..<lanes {
                    let j = s + lane * stripeCount
                    if j < qlen {
                        let sc = mat[k * m + Int(query[j])]
                        vec[lane] = UInt8(clamping: Int(sc) + Int(bias))
                    }
                }
                profile[k][s] = vec
            }
        }

        // DP vectors
        var H = [SIMD16<UInt8>](repeating: zeroVec, count: stripeCount)
        var E = [SIMD16<UInt8>](repeating: zeroVec, count: stripeCount)

        let biasVec = SIMD16<UInt8>(repeating: bias)
        let gapO = SIMD16<UInt8>(repeating: UInt8(scoring.gapOpenPenalty))
        let gapE = SIMD16<UInt8>(repeating: UInt8(scoring.gapExtendPenalty))

        var maxScore: UInt8 = UInt8(clamping: max(0, Int(h0)))
        var maxI: Int32 = -1
        var maxJ: Int32 = -1
        var overflow = false

        for i in 0..<tlen {
            let targetBase = Int(target[i])
            let prof = profile[targetBase]

            var f = SIMD16<UInt8>(repeating: 0)

            for s in 0..<stripeCount {
                // H_new = H_prev_diag + profile[target[i]][j] - bias
                var hNew = H[s] &+ prof[s]
                // Saturating subtraction of bias: if hNew < biasVec, clamp to 0
                hNew = pointwiseMax(hNew, biasVec) &- biasVec

                // E = max(H - gapO, E) - gapE (saturating)
                // Saturating: H[s] - gapO, clamped to 0
                let hSatGapO = pointwiseMax(H[s], gapO) &- gapO
                E[s] = pointwiseMax(hSatGapO, E[s])
                E[s] = pointwiseMax(E[s], gapE) &- gapE

                // F = max(H - gapO, F) - gapE (saturating)
                let hNewSatGapO = pointwiseMax(hNew, gapO) &- gapO
                f = pointwiseMax(hNewSatGapO, f)
                f = pointwiseMax(f, gapE) &- gapE

                // H = max(H_new, E, F)
                hNew = pointwiseMax(hNew, pointwiseMax(E[s], f))

                H[s] = hNew

                // Check overflow (if max approaches 255)
                let localMax = hNew.max()
                if localMax > 250 {
                    overflow = true
                }

                if localMax > maxScore {
                    maxScore = localMax
                    maxI = Int32(i)
                    // Find which lane has the max
                    for lane in 0..<lanes {
                        if hNew[lane] == localMax {
                            maxJ = Int32(s + lane * stripeCount)
                            break
                        }
                    }
                }
            }

            if overflow { return nil }
        }

        return SWResult(
            score: Int32(maxScore),
            queryEnd: maxJ + 1,
            targetEnd: maxI + 1,
            globalTargetEnd: -1,
            globalScore: -1,
            maxOff: 0
        )
    }

    @inline(__always)
    private static func pointwiseMax(_ a: SIMD16<UInt8>, _ b: SIMD16<UInt8>) -> SIMD16<UInt8> {
        a.replacing(with: b, where: b .> a)
    }
}
