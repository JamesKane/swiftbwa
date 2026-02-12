import BWACore

/// Result of a local Smith-Waterman alignment.
public struct LocalSWResult: Sendable {
    public var score: Int32
    public var queryBegin: Int32
    public var queryEnd: Int32      // inclusive
    public var targetBegin: Int32
    public var targetEnd: Int32     // inclusive

    public init(score: Int32, queryBegin: Int32, queryEnd: Int32,
                targetBegin: Int32, targetEnd: Int32) {
        self.score = score
        self.queryBegin = queryBegin
        self.queryEnd = queryEnd
        self.targetBegin = targetBegin
        self.targetEnd = targetEnd
    }
}

/// Full local Smith-Waterman with start position recovery.
///
/// Equivalent to bwa-mem2's `ksw_align2()` with `KSW_XSTART`.
/// Uses a two-pass approach:
/// 1. Forward pass: standard local SW to find best score and end positions.
/// 2. Reverse pass: run SW on reversed prefixes to recover start positions.
///
/// The inner loop uses SIMD16<UInt8> (128-bit NEON) via Farrar's striped algorithm,
/// processing 16 query positions in parallel. Falls back to scalar if scores overflow.
public struct LocalSWAligner: Sendable {

    /// Perform local Smith-Waterman alignment.
    public static func align(
        query: [UInt8], target: [UInt8], scoring: ScoringParameters,
        scoringMatrix: [Int8]? = nil
    ) -> LocalSWResult? {
        let qLen = query.count
        let tLen = target.count
        guard qLen > 0 && tLen > 0 else { return nil }

        let mat = scoringMatrix ?? scoring.scoringMatrix()
        let gapOE = scoring.gapOpenPenalty &+ scoring.gapExtendPenalty
        let gapE = scoring.gapExtendPenalty
        let bias = UInt8(scoring.mismatchPenalty)

        // Try SIMD path first (works for scores fitting in UInt8)
        let maxPossibleScore = Int32(qLen) &* scoring.matchScore
        if maxPossibleScore <= 250 {
            let simdResult = query.withUnsafeBufferPointer { qBuf in
                target.withUnsafeBufferPointer { tBuf in
                    mat.withUnsafeBufferPointer { matBuf in
                        simdLocalSW(
                            query: qBuf.baseAddress!, qLen: qLen,
                            target: tBuf.baseAddress!, tLen: tLen,
                            mat: matBuf.baseAddress!, gapOE: gapOE, gapE: gapE,
                            bias: bias
                        )
                    }
                }
            }
            if let result = simdResult { return result }
        }

        // Scalar fallback
        return scalarAlign(
            query: query, target: target,
            mat: mat, gapOE: gapOE, gapE: gapE
        )
    }

    // MARK: - SIMD16<UInt8> Farrar Striped Local SW

    /// SIMD16 local SW with two-pass start position recovery.
    /// Returns nil if scores overflow 8-bit.
    private static func simdLocalSW(
        query: UnsafePointer<UInt8>, qLen: Int,
        target: UnsafePointer<UInt8>, tLen: Int,
        mat: UnsafePointer<Int8>, gapOE: Int32, gapE: Int32,
        bias: UInt8
    ) -> LocalSWResult? {
        // Forward pass: find best score and end positions
        guard let fwd = simdSWPass(
            query: query, qLen: qLen, target: target, tLen: tLen,
            mat: mat, gapOE: gapOE, gapE: gapE, bias: bias
        ) else { return nil }

        guard fwd.score > 0 else { return nil }

        // Reverse pass: find start positions
        let revQLen = Int(fwd.qEnd) &+ 1
        let revTLen = Int(fwd.tEnd) &+ 1

        let revQ = UnsafeMutablePointer<UInt8>.allocate(capacity: revQLen)
        let revT = UnsafeMutablePointer<UInt8>.allocate(capacity: revTLen)
        defer { revQ.deallocate(); revT.deallocate() }

        for k in 0..<revQLen { revQ[k] = query[revQLen &- 1 &- k] }
        for k in 0..<revTLen { revT[k] = target[revTLen &- 1 &- k] }

        guard let rev = simdSWPass(
            query: revQ, qLen: revQLen, target: revT, tLen: revTLen,
            mat: mat, gapOE: gapOE, gapE: gapE, bias: bias
        ) else { return nil }

        return LocalSWResult(
            score: fwd.score,
            queryBegin: fwd.qEnd &- rev.qEnd,
            queryEnd: fwd.qEnd,
            targetBegin: fwd.tEnd &- rev.tEnd,
            targetEnd: fwd.tEnd
        )
    }

    /// Single-direction SIMD SW pass using Farrar's striped algorithm.
    /// Returns nil if 8-bit overflow detected.
    private static func simdSWPass(
        query: UnsafePointer<UInt8>, qLen: Int,
        target: UnsafePointer<UInt8>, tLen: Int,
        mat: UnsafePointer<Int8>, gapOE: Int32, gapE: Int32,
        bias: UInt8
    ) -> (score: Int32, qEnd: Int32, tEnd: Int32)? {
        let lanes = 16
        let sn = (qLen &+ lanes &- 1) / lanes

        // Allocate profile + DP arrays as raw pointers
        let profile = UnsafeMutablePointer<SIMD16<UInt8>>.allocate(capacity: 5 &* sn)
        let hArr = UnsafeMutablePointer<SIMD16<UInt8>>.allocate(capacity: sn)
        let eArr = UnsafeMutablePointer<SIMD16<UInt8>>.allocate(capacity: sn)
        defer {
            profile.deallocate()
            hArr.deallocate()
            eArr.deallocate()
        }

        // Build striped query profile
        let zeroVec = SIMD16<UInt8>(repeating: 0)
        for k in 0..<5 {
            let offset = k &* sn
            for s in 0..<sn {
                var vec = zeroVec
                for lane in 0..<lanes {
                    let j = s &+ lane &* sn
                    if j < qLen {
                        let sc = mat[k &* 5 &+ Int(query[j])]
                        vec[lane] = UInt8(clamping: Int(sc) &+ Int(bias))
                    }
                }
                profile[offset &+ s] = vec
            }
        }

        // Initialize H, E to zero
        for s in 0..<sn {
            hArr[s] = zeroVec
            eArr[s] = zeroVec
        }

        let biasVec = SIMD16<UInt8>(repeating: bias)
        let gapOEVec = SIMD16<UInt8>(repeating: UInt8(gapOE))
        let gapEVec = SIMD16<UInt8>(repeating: UInt8(gapE))

        var maxScore: UInt8 = 0
        var maxI: Int32 = -1
        var maxJ: Int32 = -1

        for i in 0..<tLen {
            let tBase = Int(min(target[i], 4))
            let prof = profile + tBase &* sn

            // Diagonal: shifted H from previous row's last stripe
            var vH = shiftRight(hArr[sn &- 1])
            var f = zeroVec

            for s in 0..<sn {
                // Diagonal + profile - bias (saturating clamp to 0)
                var hNew = vH &+ prof[s]
                hNew = pmax(hNew, biasVec) &- biasVec

                let hPrev = hArr[s]
                vH = hPrev

                // E: gap in query (vertical / insertion)
                let eOpen = pmax(hPrev, gapOEVec) &- gapOEVec
                var eVal = pmax(eOpen, eArr[s])
                eVal = pmax(eVal, gapEVec) &- gapEVec
                eArr[s] = eVal

                hNew = pmax(hNew, eVal)
                hNew = pmax(hNew, f)
                hArr[s] = hNew

                // F: gap in target (horizontal / deletion)
                let fOpen = pmax(hArr[s], gapOEVec) &- gapOEVec
                f = pmax(fOpen, f)
                f = pmax(f, gapEVec) &- gapEVec
            }

            // Lazy-F correction across stripe boundaries
            var corrF = shiftRight(f)
            for s in 0..<sn {
                let newH = pmax(hArr[s], corrF)
                if newH == hArr[s] { break }
                hArr[s] = newH
                let fOpen = pmax(hArr[s], gapOEVec) &- gapOEVec
                corrF = pmax(fOpen, corrF)
                corrF = pmax(corrF, gapEVec) &- gapEVec
            }

            // Track max score and position
            for s in 0..<sn {
                let localMax = hArr[s].max()
                if localMax > 250 { return nil }  // overflow
                if localMax > maxScore {
                    maxScore = localMax
                    maxI = Int32(i)
                    for lane in 0..<lanes {
                        if hArr[s][lane] == localMax {
                            let j = s &+ lane &* sn
                            if j < qLen {
                                maxJ = Int32(j)
                                break
                            }
                        }
                    }
                }
            }
        }

        return (Int32(maxScore), maxJ, maxI)
    }

    @inline(__always)
    private static func pmax(
        _ a: SIMD16<UInt8>, _ b: SIMD16<UInt8>
    ) -> SIMD16<UInt8> {
        a.replacing(with: b, where: b .> a)
    }

    @inline(__always)
    private static func shiftRight(_ v: SIMD16<UInt8>) -> SIMD16<UInt8> {
        var r = SIMD16<UInt8>(repeating: 0)
        r[1] = v[0]; r[2] = v[1]; r[3] = v[2]; r[4] = v[3]
        r[5] = v[4]; r[6] = v[5]; r[7] = v[6]; r[8] = v[7]
        r[9] = v[8]; r[10] = v[9]; r[11] = v[10]; r[12] = v[11]
        r[13] = v[12]; r[14] = v[13]; r[15] = v[14]
        return r
    }

    // MARK: - Scalar Fallback

    private static func scalarAlign(
        query: [UInt8], target: [UInt8],
        mat: [Int8], gapOE: Int32, gapE: Int32
    ) -> LocalSWResult? {
        let qLen = query.count
        let tLen = target.count

        let forwardResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                mat.withUnsafeBufferPointer { matBuf in
                    scalarSWPass(
                        query: qBuf.baseAddress!, qLen: qLen,
                        target: tBuf.baseAddress!, tLen: tLen,
                        mat: matBuf.baseAddress!, gapOE: gapOE, gapE: gapE
                    )
                }
            }
        }

        guard forwardResult.score > 0 else { return nil }

        let bestQEnd = forwardResult.qEnd
        let bestTEnd = forwardResult.tEnd
        let revQLen = Int(bestQEnd) &+ 1
        let revTLen = Int(bestTEnd) &+ 1

        let revQ = UnsafeMutablePointer<UInt8>.allocate(capacity: revQLen)
        let revT = UnsafeMutablePointer<UInt8>.allocate(capacity: revTLen)
        defer { revQ.deallocate(); revT.deallocate() }
        query.withUnsafeBufferPointer { qBuf in
            for k in 0..<revQLen { revQ[k] = qBuf[revQLen &- 1 &- k] }
        }
        target.withUnsafeBufferPointer { tBuf in
            for k in 0..<revTLen { revT[k] = tBuf[revTLen &- 1 &- k] }
        }

        let reverseResult = mat.withUnsafeBufferPointer { matBuf in
            scalarSWPass(
                query: revQ, qLen: revQLen, target: revT, tLen: revTLen,
                mat: matBuf.baseAddress!, gapOE: gapOE, gapE: gapE
            )
        }

        return LocalSWResult(
            score: forwardResult.score,
            queryBegin: bestQEnd &- reverseResult.qEnd,
            queryEnd: bestQEnd,
            targetBegin: bestTEnd &- reverseResult.tEnd,
            targetEnd: bestTEnd
        )
    }

    private static func scalarSWPass(
        query: UnsafePointer<UInt8>, qLen: Int,
        target: UnsafePointer<UInt8>, tLen: Int,
        mat: UnsafePointer<Int8>, gapOE: Int32, gapE: Int32
    ) -> (score: Int32, qEnd: Int32, tEnd: Int32) {
        let h = UnsafeMutablePointer<Int32>.allocate(capacity: tLen &+ 1)
        let e = UnsafeMutablePointer<Int32>.allocate(capacity: tLen &+ 1)
        defer { h.deallocate(); e.deallocate() }
        h.initialize(repeating: 0, count: tLen &+ 1)
        e.initialize(repeating: 0, count: tLen &+ 1)

        var bestScore: Int32 = 0
        var bestQEnd: Int32 = -1
        var bestTEnd: Int32 = -1

        for i in 0..<qLen {
            let qBase = Int(min(query[i], 4))
            let matRow = mat + qBase &* 5
            var f: Int32 = 0
            var hPrev: Int32 = 0

            for j in 0..<tLen {
                let j1 = j &+ 1
                let tBase = Int(min(target[j], 4))
                var hCur = hPrev &+ Int32(matRow[tBase])

                let fOpen = h[j1] &- gapOE
                let fExt = f &- gapE
                f = fOpen > fExt ? fOpen : fExt

                let eOpen = h[j] &- gapOE
                let eExt = e[j1] &- gapE
                let eVal = eOpen > eExt ? eOpen : eExt
                e[j1] = eVal

                hCur = hCur > f ? hCur : f
                hCur = hCur > eVal ? hCur : eVal
                hCur = hCur > 0 ? hCur : 0

                hPrev = h[j1]
                h[j1] = hCur

                if hCur > bestScore {
                    bestScore = hCur
                    bestQEnd = Int32(i)
                    bestTEnd = Int32(j)
                }
            }
        }

        return (bestScore, bestQEnd, bestTEnd)
    }
}
