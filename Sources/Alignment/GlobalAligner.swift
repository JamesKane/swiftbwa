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
    ///   - scoringMatrix: Pre-built scoring matrix (avoids re-allocation per call)
    /// - Returns: GlobalAlignResult with score, CIGAR, and edit distance
    public static func align(
        query: UnsafeBufferPointer<UInt8>,
        target: UnsafeBufferPointer<UInt8>,
        scoring: ScoringParameters,
        w: Int32,
        scoringMatrix: [Int8]? = nil
    ) -> GlobalAlignResult {
        let qlen = query.count
        let tlen = target.count

        guard qlen > 0 && tlen > 0 else {
            return GlobalAlignResult()
        }

        let m = 5  // alphabet size (A, C, G, T, N)
        let mat = scoringMatrix ?? scoring.scoringMatrix()

        // E = deletion (gap in query, consuming target): uses deletion penalties
        let oDel = scoring.gapOpenPenaltyDeletion
        let eDel = scoring.gapExtendPenaltyDeletion
        // F = insertion (gap in target, consuming query): uses insertion penalties
        let oIns = scoring.gapOpenPenalty
        let eIns = scoring.gapExtendPenalty
        // Pre-computed combined open+extend (matches bwa-mem2's oe_del, oe_ins)
        let oeDel = oDel &+ eDel
        let oeIns = oIns &+ eIns

        // Effective band width: ensure it covers the length difference
        let w = max(w, abs(Int32(qlen) &- Int32(tlen)))
        let nCol = min(qlen, Int(2 &* w &+ 1))

        // Allocate DP arrays as raw pointers — avoids Array metadata/ARC and ensures
        // truly unchecked subscript access (UnsafeMutablePointer has no _precondition).
        let negInf = Int32.min / 2  // avoid overflow on subtraction
        let qp = UnsafeMutablePointer<Int32>.allocate(capacity: m &* qlen)
        let h = UnsafeMutablePointer<Int32>.allocate(capacity: qlen &+ 1)
        let e = UnsafeMutablePointer<Int32>.allocate(capacity: qlen &+ 1)
        let z = UnsafeMutablePointer<UInt8>.allocate(capacity: tlen &* nCol)
        defer {
            qp.deallocate()
            h.deallocate()
            e.deallocate()
            z.deallocate()
        }

        // Build query profile: qp[k * qlen + j] = mat[k * m + query[j]]
        mat.withUnsafeBufferPointer { matBuf in
            let mp = matBuf.baseAddress!
            for k in 0..<m {
                for j in 0..<qlen {
                    qp[k &* qlen &+ j] = Int32(mp[k &* m &+ Int(query[j])])
                }
            }
        }

        // Initialize DP arrays
        h.initialize(repeating: negInf, count: qlen &+ 1)
        e.initialize(repeating: negInf, count: qlen &+ 1)
        z.initialize(repeating: 0, count: tlen &* nCol)

        // Initialize first row: h[j] = cost of j insertions (gap in target)
        h[0] = 0
        let initEnd = min(Int(w), qlen)
        for j in 1...initEnd {
            h[j] = 0 &- (oIns &+ eIns &* Int32(j))
        }

        // Main DP loop — uses wrapping arithmetic (&+, &-, &*) to eliminate overflow
        // traps, and ternary + max() to encourage branchless CSEL code generation.
        // Matches bwa-mem2's ksw_global2 inner loop structure.
        for i in 0..<tlen {
            let beg = max(0, Int(Int32(i) &- w))
            let end = min(qlen, Int(Int32(i) &+ w &+ 1))

            let targetBase = Int(target[i])
            let qpRow = targetBase &* qlen

            var f: Int32 = negInf
            var hDiag = h[beg]

            // Left boundary: cost of (i+1) deletions (gap in query)
            if beg == 0 {
                h[beg] = 0 &- (oeDel &+ eDel &* Int32(i))
            } else {
                h[beg] = negInf
            }

            let zi = z + (i &* nCol &- beg)

            for j in beg..<end {
                let j1 = j &+ 1

                // Score from diagonal (match/mismatch)
                var hCur = hDiag &+ qp[qpRow &+ j]

                // Save H(i-1,j) for next iteration's diagonal
                hDiag = h[j1]

                // E: deletion (gap in query, vertical move)
                // E[i,j] = max(H[i-1,j] - oDel - eDel, E[i-1,j] - eDel)
                let eOpen = hDiag &- oeDel
                let eExtend = e[j1] &- eDel
                let eFlag: UInt8 = eOpen > eExtend ? 0 : 0x01
                let eVal = eOpen > eExtend ? eOpen : eExtend
                e[j1] = eVal

                // F: insertion (gap in target, horizontal move)
                // F[i,j] = max(H[i,j-1] - oIns - eIns, F[i,j-1] - eIns)
                let fOpen = h[j] &- oeIns
                let fExtend = f &- eIns
                let fFlag: UInt8 = fOpen > fExtend ? 0 : 0x02
                let fVal = fOpen > fExtend ? fOpen : fExtend
                f = fVal

                // H = max(H_diag + score, E, F) — two-step max for correct direction
                let fromE: UInt8 = eVal > hCur ? 0x04 : 0
                hCur = eVal > hCur ? eVal : hCur
                let fromF: UInt8 = fVal > hCur ? 0x08 : 0
                hCur = fVal > hCur ? fVal : hCur

                h[j1] = hCur
                zi[j] = eFlag | fFlag | fromE | fromF
            }
        }

        let score = h[qlen]

        // Traceback from (tlen-1, qlen-1)
        let cigar = traceback(z: z, nCol: nCol, qlen: qlen, tlen: tlen, w: w)

        // NM is computed by CIGARGenerator after deletion squeeze, not here.
        return GlobalAlignResult(score: score, cigar: cigar)
    }

    /// Traceback through the backpointer matrix to produce CIGAR.
    ///
    /// Follows bwa-mem2's ksw_global2 traceback: state machine with M/E/F states.
    /// - M state: check bits 2-3 for H source. If diagonal, emit M and move diag.
    ///   If from E/F, switch state without emitting or moving.
    /// - E state (deletion): emit D, check bit 0 for extend, decrement i.
    /// - F state (insertion): emit I, check bit 1 for extend, decrement j.
    private static func traceback(
        z: UnsafePointer<UInt8>,
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

            let beg = max(0, Int(Int32(i) &- w))
            let col = j &- beg

            guard col >= 0 && col < nCol else {
                break
            }

            let d = z[i &* nCol &+ col]

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
