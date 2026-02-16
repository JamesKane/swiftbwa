import BWACore

/// Scalar banded Smith-Waterman extension.
/// Matches bwa-mem2's `scalarBandedSWA()` from bandedSWA.cpp, including:
/// - Dynamic band narrowing (shrinks beg/end when H/E are 0)
/// - W cap based on maximum possible gap (with end_bonus)
/// - Gap-cost-adjusted z-drop
/// - Diagonal restart prevention (M = M? M + q[j] : 0)
public struct BandedSWScalar: Sendable {

    /// Perform banded Smith-Waterman extension.
    ///
    /// - Parameters:
    ///   - query: Query sequence (2-bit encoded, length qlen)
    ///   - target: Target/reference sequence (2-bit encoded, length tlen)
    ///   - scoring: Scoring parameters
    ///   - w: Band width
    ///   - h0: Initial H score (score accumulated before this extension)
    ///   - endBonus: Bonus for reaching query end (pen_clip5 or pen_clip3)
    /// - Returns: SWResult with score, endpoints, and max offset
    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        h0: Int32,
        endBonus: Int32 = 0
    ) -> SWResult {
        let qlen = query.count
        let tlen = target.count
        let m = 5  // alphabet size (A,C,G,T,N)

        guard qlen > 0 && tlen > 0 else {
            return SWResult()
        }

        let mat = scoring.scoringMatrix()
        let oDel = scoring.gapOpenPenaltyDeletion
        let eDel = scoring.gapExtendPenaltyDeletion
        let oIns = scoring.gapOpenPenalty
        let eIns = scoring.gapExtendPenalty
        let zDrop = scoring.zDrop

        // Build query profile: qp[k * qlen + j] = mat[k * m + query[j]]
        var qp = [Int32](repeating: 0, count: m * qlen)
        for k in 0..<m {
            for j in 0..<qlen {
                qp[k * qlen + j] = Int32(mat[k * m + Int(query[j])])
            }
        }

        // Cap w based on maximum possible gap (bwa-mem2 bandedSWA.cpp:151-160)
        var w = w
        var maxMatScore: Int32 = 0
        for i in 0..<(m * m) {
            let s = Int32(mat[i])
            if s > maxMatScore { maxMatScore = s }
        }
        let maxIns = Int32(Double(Int32(qlen) * maxMatScore + endBonus - oIns) / Double(eIns) + 1.0)
        w = min(w, max(1, maxIns))
        let maxDel = Int32(Double(Int32(qlen) * maxMatScore + endBonus - oDel) / Double(eDel) + 1.0)
        w = min(w, max(1, maxDel))

        // DP arrays: h = H scores, e = E scores
        var h = [Int32](repeating: 0, count: qlen + 1)
        var e = [Int32](repeating: 0, count: qlen + 1)

        // Initialize first row (bwa-mem2 bandedSWA.cpp:147-149)
        h[0] = h0
        if h0 > oIns + eIns {
            h[1] = h0 - oIns - eIns
            var j = 2
            while j <= qlen && h[j - 1] > eIns {
                h[j] = h[j - 1] - eIns
                j += 1
            }
        }

        var maxScore: Int32 = h0
        var maxI: Int32 = -1
        var maxJ: Int32 = -1
        var gScore: Int32 = -1
        var gTle: Int32 = -1
        var maxOff: Int32 = 0

        var beg: Int = 0
        var end: Int = qlen

        // Main DP loop over target positions
        for i in 0..<tlen {
            let targetBase = Int(target[i])
            let qpRow = targetBase * qlen

            // Apply band constraints
            if beg < i - Int(w) { beg = i - Int(w) }
            if end > i + Int(w) + 1 { end = i + Int(w) + 1 }
            if end > qlen { end = qlen }

            var f: Int32 = 0
            // Compute first column (bwa-mem2 bandedSWA.cpp:175-178)
            var h1: Int32
            if beg == 0 {
                h1 = h0 - (oDel + eDel * Int32(i + 1))
                if h1 < 0 { h1 = 0 }
            } else {
                h1 = 0
            }

            var rowMax: Int32 = 0
            var rowMaxJ: Int32 = -1

            for j in beg..<end {
                // At loop start: h[j] = H(i-1,j-1), e[j] = E(i,j), h1 = H(i,j-1)
                let M = h[j]
                let eVal = e[j]
                h[j] = h1  // Store H(i,j-1) for next row's diagonal

                // Diagonal: M = M? M + q[j] : 0 (bwa-mem2 bandedSWA.cpp:188)
                let diag = M != 0 ? M + qp[qpRow + j] : Int32(0)

                // H(i,j) = max(M, E, F)
                var hVal = diag > eVal ? diag : eVal
                hVal = hVal > f ? hVal : f
                h1 = hVal

                // Track row max — bwa-mem2 uses >= (last position wins on ties)
                rowMaxJ = rowMax > hVal ? rowMaxJ : Int32(j)
                rowMax = rowMax > hVal ? rowMax : hVal

                // E(i+1,j) = max(M - oe_del, E(i,j) - e_del) (bwa-mem2 bandedSWA.cpp:194-198)
                var t = diag - (oDel + eDel)
                if t < 0 { t = 0 }
                var newE = eVal - eDel
                if newE < t { newE = t }
                e[j] = newE

                // F(i,j+1) = max(M - oe_ins, F(i,j) - e_ins) (bwa-mem2 bandedSWA.cpp:199-202)
                t = diag - (oIns + eIns)
                if t < 0 { t = 0 }
                var newF = f - eIns
                if newF < t { newF = t }
                f = newF
            }

            h[end] = h1
            e[end] = 0

            // Check for global score (reaching end of query) — ties go to last
            if end == qlen {
                gTle = gScore > h1 ? gTle : Int32(i)
                gScore = gScore > h1 ? gScore : h1
            }

            // Early termination when row max is 0 (bwa-mem2 bandedSWA.cpp:210)
            if rowMax == 0 { break }

            // Track overall maximum and z-drop (bwa-mem2 bandedSWA.cpp:211-220)
            if rowMax > maxScore {
                maxScore = rowMax
                maxI = Int32(i)
                maxJ = rowMaxJ
                let off = abs(rowMaxJ - Int32(i))
                if off > maxOff { maxOff = off }
            } else if zDrop > 0 {
                let deltaI = Int32(i) - maxI
                let deltaJ = rowMaxJ - maxJ
                if deltaI > deltaJ {
                    if maxScore - rowMax - (deltaI - deltaJ) * eDel > zDrop { break }
                } else {
                    if maxScore - rowMax - (deltaJ - deltaI) * eIns > zDrop { break }
                }
            }

            // Dynamic band narrowing (bwa-mem2 bandedSWA.cpp:222-225)
            var newBeg = beg
            while newBeg < end && h[newBeg] == 0 && e[newBeg] == 0 { newBeg += 1 }
            beg = newBeg
            var newEnd = end
            while newEnd >= beg && h[newEnd] == 0 && e[newEnd] == 0 { newEnd -= 1 }
            end = min(newEnd + 2, qlen)
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
}
