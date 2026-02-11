import BWACore

/// Banded global (Needleman-Wunsch) alignment with backpointer traceback.
///
/// Reimplements `ksw_global2()` from bwa-mem2's `ksw.cpp`.
/// Unlike `BandedSWScalar` (local Smith-Waterman), this performs global alignment
/// over the full query and target, producing a CIGAR string via backpointer matrix.
public struct GlobalAligner: Sendable {

    /// Perform banded global alignment with CIGAR traceback.
    ///
    /// - Parameters:
    ///   - query: 2-bit encoded query segment (qb..qe from MemAlnReg)
    ///   - target: 2-bit encoded target/reference segment (rb..re from MemAlnReg)
    ///   - scoring: Scoring parameters
    ///   - w: Band width
    /// - Returns: GlobalAlignResult with score, CIGAR, and edit distance
    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32
    ) -> GlobalAlignResult {
        let qlen = query.count
        let tlen = target.count

        guard qlen > 0 && tlen > 0 else {
            return GlobalAlignResult()
        }

        let m = 5  // alphabet size (A, C, G, T, N)
        let mat = scoring.scoringMatrix()

        // E = deletion (gap in query, consuming target): uses deletion penalties
        let oDel = scoring.gapOpenPenaltyDeletion
        let eDel = scoring.gapExtendPenaltyDeletion
        // F = insertion (gap in target, consuming query): uses insertion penalties
        let oIns = scoring.gapOpenPenalty
        let eIns = scoring.gapExtendPenalty

        // Effective band width: ensure it covers the length difference
        let w = max(w, abs(Int32(qlen) - Int32(tlen)))
        let nCol = min(qlen, Int(2 * w + 1))

        // Build query profile: qp[k * qlen + j] = mat[k * m + query[j]]
        var qp = [Int32](repeating: 0, count: m * qlen)
        for k in 0..<m {
            for j in 0..<qlen {
                qp[k * qlen + j] = Int32(mat[k * m + Int(query[j])])
            }
        }

        // DP arrays
        let negInf = Int32.min / 2  // avoid overflow on subtraction
        var h = [Int32](repeating: negInf, count: qlen + 1)
        var e = [Int32](repeating: negInf, count: qlen + 1)

        // Backpointer matrix: z[i * nCol + (j - beg)]
        // Bit 0: E source (0=opened from H, 1=extended from E)
        // Bit 1: F source (0=opened from H, 1=extended from F)
        // Bit 2: H came from E (deletion)
        // Bit 3: H came from F (insertion)
        // If neither bit 2 nor 3: H came from diagonal
        var z = [UInt8](repeating: 0, count: tlen * nCol)

        // Initialize first row: h[j] = cost of j insertions (gap in target)
        h[0] = 0
        let initEnd = min(Int(w), qlen)
        for j in 1...initEnd {
            h[j] = -(oIns + eIns * Int32(j))
            e[j] = negInf
        }

        // Main DP loop over target positions
        for i in 0..<tlen {
            let beg = max(0, Int(Int32(i) - w))
            let end = min(qlen, Int(Int32(i) + w + 1))

            let targetBase = Int(target[i])
            let qpRow = targetBase * qlen

            var f: Int32 = negInf
            var hDiag = h[beg]

            // Left boundary: cost of (i+1) deletions (gap in query)
            if beg == 0 {
                h[beg] = -(oDel + eDel * Int32(i + 1))
            } else {
                h[beg] = negInf
            }

            for j in beg..<end {
                let zIdx = i * nCol + (j - beg)
                var d: UInt8 = 0

                // Score from diagonal (match/mismatch)
                var hCur = hDiag + qp[qpRow + j]

                hDiag = h[j + 1]

                // E: deletion (gap in query, vertical move)
                // E[i,j] = max(H[i-1,j] - oDel - eDel, E[i-1,j] - eDel)
                let eOpen = h[j + 1] - oDel - eDel
                var eVal = e[j + 1] - eDel
                if eOpen > eVal {
                    eVal = eOpen
                    // bit 0 stays 0: E opened from H
                } else {
                    d |= 0x01  // E extends
                }
                e[j + 1] = eVal

                // F: insertion (gap in target, horizontal move)
                // F[i,j] = max(H[i,j-1] - oIns - eIns, F[i,j-1] - eIns)
                let fOpen = h[j] - oIns - eIns
                var fVal = f - eIns
                if fOpen > fVal {
                    fVal = fOpen
                    // bit 1 stays 0: F opened from H
                } else {
                    d |= 0x02  // F extends
                }
                f = fVal

                // H = max(H_diag + score, E, F)
                if eVal > hCur {
                    hCur = eVal
                    d |= 0x04  // H from E (deletion)
                }
                if fVal > hCur {
                    hCur = fVal
                    d |= 0x08  // H from F (insertion)
                }

                h[j + 1] = hCur
                z[zIdx] = d
            }
        }

        let score = h[qlen]

        // Traceback from (tlen-1, qlen-1)
        let cigar = traceback(z: z, nCol: nCol, qlen: qlen, tlen: tlen, w: w)

        // Compute NM from CIGAR + sequences
        let nm = computeNM(cigar: cigar, query: query, target: target)

        return GlobalAlignResult(score: score, cigar: cigar, nm: nm)
    }

    /// Traceback through the backpointer matrix to produce CIGAR.
    ///
    /// Follows bwa-mem2's ksw_global2 traceback: state machine with M/E/F states.
    /// - M state: check bits 2-3 for H source. If diagonal, emit M and move diag.
    ///   If from E/F, switch state without emitting or moving.
    /// - E state (deletion): emit D, check bit 0 for extend, decrement i.
    /// - F state (insertion): emit I, check bit 1 for extend, decrement j.
    private static func traceback(
        z: [UInt8],
        nCol: Int,
        qlen: Int,
        tlen: Int,
        w: Int32
    ) -> [UInt32] {
        var builder = CIGARBuilder()
        var i = tlen - 1  // target position
        var j = qlen - 1  // query position
        var state = 0      // 0=M, 1=E(deletion), 2=F(insertion)

        while i >= 0 || j >= 0 {
            if i < 0 {
                // Only query bases remain -> insertions
                builder.append(.insertion)
                j -= 1
                continue
            }
            if j < 0 {
                // Only target bases remain -> deletions
                builder.append(.deletion)
                i -= 1
                continue
            }

            let beg = max(0, Int(Int32(i) - w))
            let col = j - beg

            guard col >= 0 && col < nCol else {
                break
            }

            let d = z[i * nCol + col]

            switch state {
            case 0: // M state
                let hDir = Int(d >> 2) & 0x03
                if hDir == 0 {
                    // Diagonal: emit M, move both
                    builder.append(.match)
                    i -= 1
                    j -= 1
                } else if hDir == 1 {
                    // H came from E -> switch to deletion state (no emit, no move)
                    state = 1
                } else {
                    // H came from F (hDir == 2 or 3) -> switch to insertion state
                    state = 2
                }
            case 1: // E state (deletion: gap in query, consume target)
                builder.append(.deletion)
                // Check extend bit BEFORE moving
                if (d & 0x01) == 0 {
                    state = 0  // E was opened from H, return to M
                }
                i -= 1
            case 2: // F state (insertion: gap in target, consume query)
                builder.append(.insertion)
                // Check extend bit BEFORE moving
                if (d & 0x02) == 0 {
                    state = 0  // F was opened from H, return to M
                }
                j -= 1
            default:
                break
            }
        }

        // Traceback produces CIGAR in reverse order
        builder.reverse()
        return builder.build()
    }

    /// Compute edit distance (NM) by walking CIGAR and comparing bases.
    static func computeNM(
        cigar: [UInt32],
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>
    ) -> Int32 {
        var nm: Int32 = 0
        var qi = 0
        var ti = 0

        for c in cigar {
            let op = c & 0xF
            let len = Int(c >> 4)

            switch op {
            case CIGAROp.match.rawValue:
                for _ in 0..<len {
                    if qi < query.count && ti < target.count {
                        if query[qi] != target[ti] {
                            nm += 1
                        }
                        qi += 1
                        ti += 1
                    }
                }
            case CIGAROp.insertion.rawValue:
                nm += Int32(len)
                qi += len
            case CIGAROp.deletion.rawValue:
                nm += Int32(len)
                ti += len
            case CIGAROp.softClip.rawValue:
                qi += len
            default:
                break
            }
        }

        return nm
    }
}
