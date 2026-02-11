import BWACore

/// 16-bit SIMD banded Smith-Waterman using SIMD8<Int16> (128-bit NEON).
/// Used when scores don't fit in 8-bit range.
public struct BandedSW16: Sendable {

    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        h0: Int32
    ) -> SWResult {
        let qlen = query.count
        let tlen = target.count

        guard qlen > 0 && tlen > 0 else { return SWResult() }

        let mat = scoring.scoringMatrix()
        let m = 5
        let lanes = 8  // SIMD8<Int16> width
        let stripeCount = (qlen + lanes - 1) / lanes

        guard stripeCount > 0 else { return SWResult() }

        // Build striped query profile
        var profile = [[SIMD8<Int16>]](repeating: [SIMD8<Int16>](repeating: .zero, count: stripeCount), count: m)

        for k in 0..<m {
            for s in 0..<stripeCount {
                var vec = SIMD8<Int16>(repeating: 0)
                for lane in 0..<lanes {
                    let j = s + lane * stripeCount
                    if j < qlen {
                        vec[lane] = Int16(mat[k * m + Int(query[j])])
                    }
                }
                profile[k][s] = vec
            }
        }

        var H = [SIMD8<Int16>](repeating: .zero, count: stripeCount)
        var E = [SIMD8<Int16>](repeating: .zero, count: stripeCount)

        let gapO = SIMD8<Int16>(repeating: Int16(scoring.gapOpenPenalty))
        let gapE = SIMD8<Int16>(repeating: Int16(scoring.gapExtendPenalty))
        let zero = SIMD8<Int16>(repeating: 0)

        var maxScore: Int16 = Int16(clamping: Int(h0))
        var maxI: Int32 = -1
        var maxJ: Int32 = -1

        for i in 0..<tlen {
            let targetBase = Int(target[i])
            let prof = profile[targetBase]

            var f = SIMD8<Int16>(repeating: 0)

            for s in 0..<stripeCount {
                var hNew = H[s] &+ prof[s]
                hNew = hNew.replacing(with: zero, where: hNew .< zero)

                // E
                let hMinusGap = H[s] &- gapO
                E[s] = pointwiseMax(hMinusGap, E[s]) &- gapE
                E[s] = E[s].replacing(with: zero, where: E[s] .< zero)

                // F
                let hForF = hNew &- gapO
                f = pointwiseMax(hForF, f) &- gapE
                f = f.replacing(with: zero, where: f .< zero)

                hNew = pointwiseMax(hNew, pointwiseMax(E[s], f))
                H[s] = hNew

                let localMax = hNew.max()
                if localMax > maxScore {
                    maxScore = localMax
                    maxI = Int32(i)
                    for lane in 0..<lanes {
                        if hNew[lane] == localMax {
                            maxJ = Int32(s + lane * stripeCount)
                            break
                        }
                    }
                }
            }
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
    private static func pointwiseMax(_ a: SIMD8<Int16>, _ b: SIMD8<Int16>) -> SIMD8<Int16> {
        a.replacing(with: b, where: b .> a)
    }
}
