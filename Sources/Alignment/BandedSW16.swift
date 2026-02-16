import BWACore

/// 16-bit SIMD banded Smith-Waterman using SIMD8<Int16> (128-bit NEON).
/// Implements Farrar's striped algorithm with z-dropoff and lazy-F correction.
/// Matches bwa-mem2's ksw_extl() behavior.
public struct BandedSW16: Sendable {

    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        h0: Int32,
        scoringMatrix: [Int8]? = nil
    ) -> SWResult {
        let qlen = query.count
        let tlen = target.count

        guard qlen > 0 && tlen > 0 else { return SWResult() }

        let mat = scoringMatrix ?? scoring.scoringMatrix()
        let m = 5
        let lanes = 8  // SIMD8<Int16> width
        let stripeCount = (qlen + lanes - 1) / lanes

        guard stripeCount > 0 else { return SWResult() }

        let oIns = scoring.gapOpenPenalty
        let eIns = scoring.gapExtendPenalty
        let oDel = scoring.gapOpenPenaltyDeletion
        let eDel = scoring.gapExtendPenaltyDeletion
        let zDrop = scoring.zDrop

        let zero = SIMD8<Int16>(repeating: 0)

        // Single allocation for all DP arrays:
        // profile: m * stripeCount, H: stripeCount, E: stripeCount
        let totalVecs = m * stripeCount + stripeCount + stripeCount
        let mem = UnsafeMutablePointer<SIMD8<Int16>>.allocate(capacity: totalVecs)
        mem.initialize(repeating: zero, count: totalVecs)
        defer { mem.deallocate() }

        let profile = mem                                        // [0 ..< m*stripeCount]
        let H = mem.advanced(by: m * stripeCount)                // [m*stripeCount ..< m*stripeCount + stripeCount]
        let E = mem.advanced(by: m * stripeCount + stripeCount)  // [... + stripeCount ..< totalVecs]

        // Build striped query profile
        for k in 0..<m {
            let profK = profile.advanced(by: k * stripeCount)
            for s in 0..<stripeCount {
                var vec = SIMD8<Int16>(repeating: 0)
                for lane in 0..<lanes {
                    let j = s + lane * stripeCount
                    if j < qlen {
                        vec[lane] = Int16(mat[k * m + Int(query[j])])
                    }
                }
                profK[s] = vec
            }
        }

        // Initialize H: H[position p] = max(0, h0 - oIns - eIns*(p+1)) for p+1 <= w
        for s in 0..<stripeCount {
            var vec = SIMD8<Int16>(repeating: 0)
            for lane in 0..<lanes {
                let j = s + lane * stripeCount
                if j < qlen && j + 1 <= Int(w) {
                    let val = Int(h0) - Int(oIns) - Int(eIns) * (j + 1)
                    if val > 0 { vec[lane] = Int16(clamping: val) }
                }
            }
            H[s] = vec
        }

        let oInsVec = SIMD8<Int16>(repeating: Int16(oIns))
        let eInsVec = SIMD8<Int16>(repeating: Int16(eIns))
        let oDelVec = SIMD8<Int16>(repeating: Int16(oDel))
        let eDelVec = SIMD8<Int16>(repeating: Int16(eDel))

        var maxScore: Int32 = h0
        var maxI: Int32 = -1
        var maxJ: Int32 = -1
        var gScore: Int32 = -1
        var gTle: Int32 = -1
        var maxOff: Int32 = 0

        for i in 0..<tlen {
            let targetBase = Int(target[i])
            let prof = profile.advanced(by: targetBase * stripeCount)

            // Diagonal: shifted H from previous row's last stripe
            var vH = shiftLanesRight(H[stripeCount - 1])
            if i == 0 {
                // Inject h0 as the diagonal for position 0 on the first target row
                vH[0] = Int16(clamping: max(0, Int(h0)))
            }

            var f = SIMD8<Int16>(repeating: 0)

            for s in 0..<stripeCount {
                // Diagonal contribution
                var hNew = vH &+ prof[s]
                hNew = hNew.replacing(with: zero, where: hNew .< zero)
                // Zero out diagonal contribution where vH was 0 (prevents alignment restart).
                // Matches bwa-mem2's ksw_extend2: M = M? M + q[j] : 0
                hNew = hNew.replacing(with: zero, where: vH .== zero)

                // Save H_prev[s] and set vH for next stripe's diagonal
                let hPrev = H[s]
                vH = hPrev

                // E: insertion from previous row
                let hMinusGapIns = hPrev &- oInsVec
                E[s] = pointwiseMax(hMinusGapIns, E[s]) &- eInsVec
                E[s] = E[s].replacing(with: zero, where: E[s] .< zero)

                // Combine diagonal with E
                hNew = pointwiseMax(hNew, E[s])

                // Combine with F (deletion, carried from previous stripe)
                hNew = pointwiseMax(hNew, f)

                // Write new H
                H[s] = hNew

                // Compute F for next stripe using final H
                let hForF = H[s] &- oDelVec
                f = pointwiseMax(hForF, f) &- eDelVec
                f = f.replacing(with: zero, where: f .< zero)
            }

            // Lazy-F correction: handle F propagation across lane group boundaries
            var corrF = shiftLanesRight(f)
            for s in 0..<stripeCount {
                let newH = pointwiseMax(H[s], corrF)
                if newH == H[s] { break }
                H[s] = newH
                let hForF = H[s] &- oDelVec
                corrF = pointwiseMax(hForF, corrF) &- eDelVec
                corrF = corrF.replacing(with: zero, where: corrF .< zero)
            }

            // Tracking pass: row max (and its column), overflow check
            var rowMax: Int32 = 0
            var rowMaxJ: Int32 = -1
            for s in 0..<stripeCount {
                let localMax = H[s].max()
                let localMax32 = Int32(localMax)
                if localMax32 > rowMax {
                    rowMax = localMax32
                    for lane in 0..<lanes {
                        if H[s][lane] == localMax {
                            let j = s + lane * stripeCount
                            if j < qlen {
                                rowMaxJ = Int32(j)
                                break
                            }
                        }
                    }
                }
            }

            // Early termination when row max is 0 (bwa-mem2 ksw.cpp:511)
            if rowMax == 0 { break }

            // Track overall maximum and z-drop
            if rowMax > maxScore {
                maxScore = rowMax
                maxI = Int32(i)
                maxJ = rowMaxJ
                let off = abs(rowMaxJ - Int32(i))
                if off > maxOff { maxOff = off }
            } else if zDrop > 0 {
                // Z-drop: gap-cost-adjusted formula (bwa-mem2 ksw.cpp:515-520)
                let deltaI = Int32(i) - maxI
                let deltaJ = rowMaxJ - maxJ
                if deltaI > deltaJ {
                    if maxScore - rowMax - (deltaI - deltaJ) * eDel > zDrop { break }
                } else {
                    if maxScore - rowMax - (deltaJ - deltaI) * eIns > zDrop { break }
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
        }

        return SWResult(
            score: maxScore,
            queryEnd: maxJ + 1,
            targetEnd: maxI + 1,
            globalTargetEnd: gTle + 1,
            globalScore: gScore,
            maxOff: maxOff
        )
    }

    @inline(__always)
    private static func pointwiseMax(_ a: SIMD8<Int16>, _ b: SIMD8<Int16>) -> SIMD8<Int16> {
        a.replacing(with: b, where: b .> a)
    }

    /// Shift SIMD lanes right by one: result[0] = 0, result[i] = source[i-1]
    @inline(__always)
    private static func shiftLanesRight(_ v: SIMD8<Int16>) -> SIMD8<Int16> {
        var r = SIMD8<Int16>(repeating: 0)
        r[1] = v[0]; r[2] = v[1]; r[3] = v[2]; r[4] = v[3]
        r[5] = v[4]; r[6] = v[5]; r[7] = v[6]
        return r
    }
}
