import BWACore

/// Scalar (non-SIMD) banded Smith-Waterman alignment.
/// Reference implementation matching bwa-mem2's `scalarBandedSWA()`.
public struct BandedSWScalar: Sendable {

    /// Perform banded Smith-Waterman extension.
    ///
    /// - Parameters:
    ///   - query: Query sequence (2-bit encoded, length qlen)
    ///   - target: Target/reference sequence (2-bit encoded, length tlen)
    ///   - scoring: Scoring parameters
    ///   - w: Band width
    ///   - h0: Initial H score (score accumulated before this extension)
    /// - Returns: SWResult with score, endpoints, and max offset
    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        h0: Int32
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

        // DP arrays: h = H scores, e = E scores (insertion)
        var h = [Int32](repeating: 0, count: qlen + 1)
        var e = [Int32](repeating: 0, count: qlen + 1)

        // Initialize
        h[0] = h0
        let initEnd = min(Int(w), qlen)
        for j in 1...initEnd {
            h[j] = h0 > oIns + eIns * Int32(j) ? h0 - oIns - eIns * Int32(j) : 0
            e[j] = 0
        }
        if initEnd + 1 <= qlen {
            for j in (initEnd + 1)...qlen {
                h[j] = 0
                e[j] = 0
            }
        }

        var maxScore: Int32 = h0
        var maxI: Int32 = -1
        var maxJ: Int32 = -1
        var maxIE: Int32 = -1  // max offset tracking
        var gScore: Int32 = -1
        var gTle: Int32 = -1

        // Main DP loop over target positions
        for i in 0..<tlen {
            let beg = max(0, Int(Int32(i) - w))
            let end = min(qlen, Int(Int32(i) + w + 1))

            let targetBase = Int(target[i])
            let qpRow = targetBase * qlen  // offset into query profile

            var f: Int32 = 0  // F score (deletion)
            var hDiag = h[beg]  // H[i-1, j-1]

            h[beg] = 0

            for j in beg..<end {
                // Score from diagonal (match/mismatch)
                var hCur = hDiag + qp[qpRow + j]
                if hCur < 0 { hCur = 0 }

                hDiag = h[j + 1]

                // E = max(H - o_ins, E) - e_ins (insertion in query)
                var eVal = e[j + 1]
                let hMinusGapIns = h[j + 1] - oIns
                eVal = max(hMinusGapIns, eVal) - eIns
                if eVal < 0 { eVal = 0 }
                e[j + 1] = eVal

                // F = max(H - o_del, F) - e_del (deletion from reference)
                let hMinusGapDel = h[j] - oDel
                f = max(hMinusGapDel, f) - eDel
                if f < 0 { f = 0 }

                // H = max(H_diag + score, E, F)
                hCur = max(hCur, max(eVal, f))
                h[j + 1] = hCur

                // Track maximum score
                if hCur > maxScore {
                    maxScore = hCur
                    maxI = Int32(i)
                    maxJ = Int32(j)
                    maxIE = Int32(i)
                }

            }

            // Check for global alignment (reaching end of query)
            if end == qlen && h[qlen] > gScore {
                gScore = h[qlen]
                gTle = Int32(i)
            }

            // Z-drop: terminate if score has dropped too far
            if maxScore - h[end] > zDrop {
                // Only break if we've gone past the max
                if Int32(i) - maxIE > Int32(w) {
                    break
                }
            }
        }

        return SWResult(
            score: maxScore,
            queryEnd: maxJ + 1,
            targetEnd: maxI + 1,
            globalTargetEnd: gTle + 1,
            globalScore: gScore
        )
    }
}
