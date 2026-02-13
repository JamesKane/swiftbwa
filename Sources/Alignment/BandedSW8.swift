import BWACore

/// 8-bit SIMD banded Smith-Waterman using SIMD16<UInt8> (128-bit NEON).
/// Implements Farrar's striped algorithm with z-dropoff and lazy-F correction.
/// Falls back to nil if scores overflow UInt8 range.
public struct BandedSW8: Sendable {

    /// Perform 8-bit SIMD banded Smith-Waterman.
    /// Returns nil if scores overflow 8-bit range (caller should use BandedSW16).
    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        h0: Int32,
        scoringMatrix: [Int8]? = nil
    ) -> SWResult? {
        let qlen = query.count
        let tlen = target.count

        guard qlen > 0 && tlen > 0 else { return SWResult() }

        // Check if scores can fit in 8-bit
        let maxPossibleScore = h0 + Int32(qlen) * scoring.matchScore
        if maxPossibleScore > 250 {
            return nil  // Overflow risk, use 16-bit
        }

        let mat = scoringMatrix ?? scoring.scoringMatrix()
        let m = 5  // alphabet size
        let lanes = 16  // SIMD16<UInt8> width
        let stripeCount = (qlen + lanes - 1) / lanes

        guard stripeCount > 0 else { return SWResult() }

        let oIns = scoring.gapOpenPenalty
        let eIns = scoring.gapExtendPenalty
        let oDel = scoring.gapOpenPenaltyDeletion
        let eDel = scoring.gapExtendPenaltyDeletion
        let zDrop = scoring.zDrop

        let bias = UInt8(scoring.mismatchPenalty)
        let zeroVec = SIMD16<UInt8>(repeating: 0)

        // Single allocation for all DP arrays:
        // profile: m * stripeCount, H: stripeCount, E: stripeCount
        let totalVecs = m * stripeCount + stripeCount + stripeCount
        let mem = UnsafeMutablePointer<SIMD16<UInt8>>.allocate(capacity: totalVecs)
        mem.initialize(repeating: zeroVec, count: totalVecs)
        defer { mem.deallocate() }

        let profile = mem                                        // [0 ..< m*stripeCount]
        let H = mem.advanced(by: m * stripeCount)                // [m*stripeCount ..< m*stripeCount + stripeCount]
        let E = mem.advanced(by: m * stripeCount + stripeCount)  // [... + stripeCount ..< totalVecs]

        // Build striped query profile (with bias added to keep values unsigned)
        for k in 0..<m {
            let profK = profile.advanced(by: k * stripeCount)
            for s in 0..<stripeCount {
                var vec = SIMD16<UInt8>(repeating: 0)
                for lane in 0..<lanes {
                    let j = s + lane * stripeCount
                    if j < qlen {
                        let sc = mat[k * m + Int(query[j])]
                        vec[lane] = UInt8(clamping: Int(sc) + Int(bias))
                    }
                }
                profK[s] = vec
            }
        }

        // Initialize H: H[position p] = max(0, h0 - oIns - eIns*(p+1)) for p+1 <= w
        for s in 0..<stripeCount {
            var vec = SIMD16<UInt8>(repeating: 0)
            for lane in 0..<lanes {
                let j = s + lane * stripeCount
                if j < qlen && j + 1 <= Int(w) {
                    let val = Int(h0) - Int(oIns) - Int(eIns) * (j + 1)
                    if val > 0 { vec[lane] = UInt8(clamping: val) }
                }
            }
            H[s] = vec
        }

        let biasVec = SIMD16<UInt8>(repeating: bias)
        let oInsVec = SIMD16<UInt8>(repeating: UInt8(oIns))
        let eInsVec = SIMD16<UInt8>(repeating: UInt8(eIns))
        let oDelVec = SIMD16<UInt8>(repeating: UInt8(oDel))
        let eDelVec = SIMD16<UInt8>(repeating: UInt8(eDel))

        var maxScore: UInt8 = UInt8(clamping: max(0, Int(h0)))
        var maxI: Int32 = -1
        var maxJ: Int32 = -1
        var maxIE: Int32 = -1
        var gScore: Int32 = -1
        var gTle: Int32 = -1
        var maxOff: Int32 = 0

        for i in 0..<tlen {
            let targetBase = Int(target[i])
            let prof = profile.advanced(by: targetBase * stripeCount)

            // Diagonal: shifted H from previous row's last stripe
            var vH = shiftLanesRight(H[stripeCount - 1])
            if i == 0 {
                vH[0] = UInt8(clamping: max(0, Int(h0)))
            }

            var f = SIMD16<UInt8>(repeating: 0)

            for s in 0..<stripeCount {
                // Diagonal + profile (includes bias), then subtract bias
                var hNew = vH &+ prof[s]
                hNew = pointwiseMax(hNew, biasVec) &- biasVec

                // Save H_prev[s] and set vH for next stripe's diagonal
                let hPrev = H[s]
                vH = hPrev

                // E: insertion from previous row (saturating unsigned subtraction)
                let hSatGapOIns = pointwiseMax(hPrev, oInsVec) &- oInsVec
                E[s] = pointwiseMax(hSatGapOIns, E[s])
                E[s] = pointwiseMax(E[s], eInsVec) &- eInsVec

                // Combine diagonal with E
                hNew = pointwiseMax(hNew, E[s])

                // Combine with F (deletion, carried from previous stripe)
                hNew = pointwiseMax(hNew, f)

                // Write new H
                H[s] = hNew

                // Compute F for next stripe using final H (saturating)
                let hSatGapODel = pointwiseMax(H[s], oDelVec) &- oDelVec
                f = pointwiseMax(hSatGapODel, f)
                f = pointwiseMax(f, eDelVec) &- eDelVec
            }

            // Lazy-F correction
            var corrF = shiftLanesRight(f)
            for s in 0..<stripeCount {
                let newH = pointwiseMax(H[s], corrF)
                if newH == H[s] { break }
                H[s] = newH
                let hSatGapODel = pointwiseMax(H[s], oDelVec) &- oDelVec
                corrF = pointwiseMax(hSatGapODel, corrF)
                corrF = pointwiseMax(corrF, eDelVec) &- eDelVec
            }

            // Tracking pass: check overflow, row max, overall max, global score
            var rowMax: UInt8 = 0
            for s in 0..<stripeCount {
                let localMax = H[s].max()
                if localMax > 250 { return nil }  // Overflow
                if localMax > rowMax { rowMax = localMax }
                if localMax > maxScore {
                    maxScore = localMax
                    maxI = Int32(i)
                    maxIE = Int32(i)
                    for lane in 0..<lanes {
                        if H[s][lane] == localMax {
                            let j = s + lane * stripeCount
                            if j < qlen {
                                maxJ = Int32(j)
                                break
                            }
                        }
                    }
                    let off = abs(maxI - maxJ)
                    if off > maxOff { maxOff = off }
                }
            }

            // Track global score (alignment consuming entire query)
            let lastJ = qlen - 1
            let lastStripe = lastJ % stripeCount
            let lastLane = lastJ / stripeCount
            let hAtEnd = Int32(H[lastStripe][lastLane])
            if hAtEnd > gScore {
                gScore = hAtEnd
                gTle = Int32(i)
            }

            // Z-dropoff
            if Int32(maxScore) - Int32(rowMax) > zDrop && Int32(i) - maxIE > w {
                break
            }
        }

        return SWResult(
            score: Int32(maxScore),
            queryEnd: maxJ + 1,
            targetEnd: maxI + 1,
            globalTargetEnd: gTle + 1,
            globalScore: gScore,
            maxOff: maxOff
        )
    }

    @inline(__always)
    private static func pointwiseMax(_ a: SIMD16<UInt8>, _ b: SIMD16<UInt8>) -> SIMD16<UInt8> {
        a.replacing(with: b, where: b .> a)
    }

    /// Shift SIMD lanes right by one: result[0] = 0, result[i] = source[i-1]
    @inline(__always)
    private static func shiftLanesRight(_ v: SIMD16<UInt8>) -> SIMD16<UInt8> {
        var r = SIMD16<UInt8>(repeating: 0)
        r[1] = v[0]; r[2] = v[1]; r[3] = v[2]; r[4] = v[3]
        r[5] = v[4]; r[6] = v[5]; r[7] = v[6]; r[8] = v[7]
        r[9] = v[8]; r[10] = v[9]; r[11] = v[10]; r[12] = v[11]
        r[13] = v[12]; r[14] = v[13]; r[15] = v[14]
        return r
    }
}
