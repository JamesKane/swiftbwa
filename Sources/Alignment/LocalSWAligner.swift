import BWACore

/// Result of a local Smith-Waterman alignment.
public struct LocalSWResult: Sendable {
    public var score: Int32
    public var queryBegin: Int32
    public var queryEnd: Int32      // inclusive
    public var targetBegin: Int32
    public var targetEnd: Int32     // inclusive
    /// Second-best alignment score from a peak outside the primary alignment zone.
    /// Used as csub for tandem repeat detection. 0 if no second peak found.
    public var score2: Int32

    public init(score: Int32, queryBegin: Int32, queryEnd: Int32,
                targetBegin: Int32, targetEnd: Int32, score2: Int32 = 0) {
        self.score = score
        self.queryBegin = queryBegin
        self.queryEnd = queryEnd
        self.targetBegin = targetBegin
        self.targetEnd = targetEnd
        self.score2 = score2
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
    ///
    /// - Parameters:
    ///   - minSubScore: Minimum score threshold for tracking second-best peaks (score2/csub).
    ///     When > 0, enables tandem repeat detection matching bwa-mem2's KSW_XSUBO.
    ///     Typically `minSeedLength * matchScore`.
    public static func align(
        query: [UInt8], target: [UInt8], scoring: ScoringParameters,
        scoringMatrix: [Int8]? = nil, minSubScore: UInt8 = 0
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
                            bias: bias, minSubScore: minSubScore
                        )
                    }
                }
            }
            if let result = simdResult { return result }
        }

        // Scalar fallback
        return scalarAlign(
            query: query, target: target,
            mat: mat, gapOE: gapOE, gapE: gapE,
            minSubScore: minSubScore
        )
    }

    // MARK: - SIMD16<UInt8> Farrar Striped Local SW

    /// SIMD16 local SW with two-pass start position recovery.
    /// Returns nil if scores overflow 8-bit.
    private static func simdLocalSW(
        query: UnsafePointer<UInt8>, qLen: Int,
        target: UnsafePointer<UInt8>, tLen: Int,
        mat: UnsafePointer<Int8>, gapOE: Int32, gapE: Int32,
        bias: UInt8, minSubScore: UInt8
    ) -> LocalSWResult? {
        // Forward pass: find best score and end positions + score2
        guard let fwd = simdSWPass(
            query: query, qLen: qLen, target: target, tLen: tLen,
            mat: mat, gapOE: gapOE, gapE: gapE, bias: bias,
            minSubScore: minSubScore
        ) else { return nil }

        guard fwd.score > 0 else { return nil }

        // Reverse pass: find start positions (no score2 needed)
        let revQLen = Int(fwd.qEnd) &+ 1
        let revTLen = Int(fwd.tEnd) &+ 1

        let revQ = UnsafeMutablePointer<UInt8>.allocate(capacity: revQLen)
        let revT = UnsafeMutablePointer<UInt8>.allocate(capacity: revTLen)
        defer { revQ.deallocate(); revT.deallocate() }

        for k in 0..<revQLen { revQ[k] = query[revQLen &- 1 &- k] }
        for k in 0..<revTLen { revT[k] = target[revTLen &- 1 &- k] }

        guard let rev = simdSWPass(
            query: revQ, qLen: revQLen, target: revT, tLen: revTLen,
            mat: mat, gapOE: gapOE, gapE: gapE, bias: bias,
            minSubScore: 0
        ) else { return nil }

        return LocalSWResult(
            score: fwd.score,
            queryBegin: fwd.qEnd &- rev.qEnd,
            queryEnd: fwd.qEnd,
            targetBegin: fwd.tEnd &- rev.tEnd,
            targetEnd: fwd.tEnd,
            score2: fwd.score2
        )
    }

    /// Single-direction SIMD SW pass using Farrar's striped algorithm.
    /// Returns nil if 8-bit overflow detected.
    ///
    /// When `minSubScore > 0`, tracks per-column peak scores to compute `score2`:
    /// the best alignment score from a peak outside the primary alignment zone.
    /// Matches bwa-mem2's KSW_XSUBO logic in ksw.cpp.
    private static func simdSWPass(
        query: UnsafePointer<UInt8>, qLen: Int,
        target: UnsafePointer<UInt8>, tLen: Int,
        mat: UnsafePointer<Int8>, gapOE: Int32, gapE: Int32,
        bias: UInt8, minSubScore: UInt8
    ) -> (score: Int32, qEnd: Int32, tEnd: Int32, score2: Int32)? {
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

        // Peak tracking for score2 (matches bwa-mem2's b[] array in ksw_u8)
        // Each entry: (score << 32) | column. Consecutive columns merged.
        var nPeaks = 0
        var mPeaks = 0
        var peaks: UnsafeMutablePointer<UInt64>? = nil
        let trackPeaks = minSubScore > 0
        defer { peaks?.deallocate() }

        // Find max base score for exclusion zone calculation
        var maxBaseScore: UInt8 = 0
        if trackPeaks {
            for k in 0..<25 {
                let v = UInt8(clamping: Int(mat[k]) &+ Int(bias))
                if v > maxBaseScore { maxBaseScore = v }
            }
            maxBaseScore = maxBaseScore > bias ? maxBaseScore &- bias : 0
        }

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

            // Compute column max and track best position
            var colMax: UInt8 = 0
            for s in 0..<sn {
                let localMax = hArr[s].max()
                if localMax > 250 { return nil }  // overflow
                if localMax > colMax { colMax = localMax }
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

            // Record peaks for score2 (bwa-mem2 ksw.cpp lines 197-205)
            if trackPeaks && colMax >= minSubScore {
                let colMaxU64 = UInt64(colMax)
                if nPeaks == 0 || Int32(peaks![nPeaks &- 1] & 0xFFFF_FFFF) &+ 1 != Int32(i) {
                    // New peak (non-consecutive column)
                    if nPeaks == mPeaks {
                        let newCap = mPeaks == 0 ? 8 : mPeaks &<< 1
                        let newBuf = UnsafeMutablePointer<UInt64>.allocate(capacity: newCap)
                        if let old = peaks {
                            newBuf.update(from: old, count: nPeaks)
                            old.deallocate()
                        }
                        peaks = newBuf
                        mPeaks = newCap
                    }
                    peaks![nPeaks] = (colMaxU64 &<< 32) | UInt64(UInt32(bitPattern: Int32(i)))
                    nPeaks &+= 1
                } else if peaks![nPeaks &- 1] >> 32 < colMaxU64 {
                    // Update last peak (consecutive column, higher score)
                    peaks![nPeaks &- 1] = (colMaxU64 &<< 32) | UInt64(UInt32(bitPattern: Int32(i)))
                }
            }
        }

        // Compute score2: best peak outside the exclusion zone around the primary hit
        var score2: Int32 = 0
        if trackPeaks, let peakBuf = peaks, maxScore > 0, maxBaseScore > 0 {
            let exclusionRadius = (Int32(maxScore) &+ Int32(maxBaseScore) &- 1) / Int32(maxBaseScore)
            let low = maxI &- exclusionRadius
            let high = maxI &+ exclusionRadius
            for p in 0..<nPeaks {
                let entry = peakBuf[p]
                let col = Int32(bitPattern: UInt32(entry & 0xFFFF_FFFF))
                let sc = Int32(entry >> 32)
                if (col < low || col > high) && sc > score2 {
                    score2 = sc
                }
            }
        }

        return (Int32(maxScore), maxJ, maxI, score2)
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
        mat: [Int8], gapOE: Int32, gapE: Int32,
        minSubScore: UInt8
    ) -> LocalSWResult? {
        let qLen = query.count
        let tLen = target.count

        let forwardResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                mat.withUnsafeBufferPointer { matBuf in
                    scalarSWPass(
                        query: qBuf.baseAddress!, qLen: qLen,
                        target: tBuf.baseAddress!, tLen: tLen,
                        mat: matBuf.baseAddress!, gapOE: gapOE, gapE: gapE,
                        minSubScore: Int32(minSubScore)
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
                mat: matBuf.baseAddress!, gapOE: gapOE, gapE: gapE,
                minSubScore: 0
            )
        }

        return LocalSWResult(
            score: forwardResult.score,
            queryBegin: bestQEnd &- reverseResult.qEnd,
            queryEnd: bestQEnd,
            targetBegin: bestTEnd &- reverseResult.tEnd,
            targetEnd: bestTEnd,
            score2: forwardResult.score2
        )
    }

    private static func scalarSWPass(
        query: UnsafePointer<UInt8>, qLen: Int,
        target: UnsafePointer<UInt8>, tLen: Int,
        mat: UnsafePointer<Int8>, gapOE: Int32, gapE: Int32,
        minSubScore: Int32
    ) -> (score: Int32, qEnd: Int32, tEnd: Int32, score2: Int32) {
        // Note: scalar path uses query as outer loop (transposed vs bwa-mem2).
        // Peak tracking uses query-row index for exclusion zone.
        let h = UnsafeMutablePointer<Int32>.allocate(capacity: tLen &+ 1)
        let e = UnsafeMutablePointer<Int32>.allocate(capacity: tLen &+ 1)
        defer { h.deallocate(); e.deallocate() }
        h.initialize(repeating: 0, count: tLen &+ 1)
        e.initialize(repeating: 0, count: tLen &+ 1)

        var bestScore: Int32 = 0
        var bestQEnd: Int32 = -1
        var bestTEnd: Int32 = -1

        // Peak tracking for score2
        let trackPeaks = minSubScore > 0
        var nPeaks = 0
        var mPeaks = 0
        var peaks: UnsafeMutablePointer<UInt64>? = nil
        var maxBaseScore: Int32 = 0
        if trackPeaks {
            for k in 0..<25 {
                let v = Int32(mat[k])
                if v > maxBaseScore { maxBaseScore = v }
            }
        }
        defer { peaks?.deallocate() }

        for i in 0..<qLen {
            let qBase = Int(min(query[i], 4))
            let matRow = mat + qBase &* 5
            var f: Int32 = 0
            var hPrev: Int32 = 0
            var rowMax: Int32 = 0

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

                if hCur > rowMax { rowMax = hCur }
                if hCur > bestScore {
                    bestScore = hCur
                    bestQEnd = Int32(i)
                    bestTEnd = Int32(j)
                }
            }

            // Record peaks (using query row as the axis, analogous to target column in SIMD path)
            if trackPeaks && rowMax >= minSubScore {
                let rmU64 = UInt64(UInt32(bitPattern: rowMax))
                let iU64 = UInt64(UInt32(bitPattern: Int32(i)))
                if nPeaks == 0 || Int32(peaks![nPeaks &- 1] & 0xFFFF_FFFF) &+ 1 != Int32(i) {
                    if nPeaks == mPeaks {
                        let newCap = mPeaks == 0 ? 8 : mPeaks &<< 1
                        let newBuf = UnsafeMutablePointer<UInt64>.allocate(capacity: newCap)
                        if let old = peaks {
                            newBuf.update(from: old, count: nPeaks)
                            old.deallocate()
                        }
                        peaks = newBuf
                        mPeaks = newCap
                    }
                    peaks![nPeaks] = (rmU64 &<< 32) | iU64
                    nPeaks &+= 1
                } else if peaks![nPeaks &- 1] >> 32 < rmU64 {
                    peaks![nPeaks &- 1] = (rmU64 &<< 32) | iU64
                }
            }
        }

        // Compute score2
        var score2: Int32 = 0
        if trackPeaks, let peakBuf = peaks, bestScore > 0, maxBaseScore > 0 {
            let exclusionRadius = (bestScore &+ maxBaseScore &- 1) / maxBaseScore
            let bestAxis = bestQEnd  // scalar uses query as outer loop
            let low = bestAxis &- exclusionRadius
            let high = bestAxis &+ exclusionRadius
            for p in 0..<nPeaks {
                let entry = peakBuf[p]
                let col = Int32(bitPattern: UInt32(entry & 0xFFFF_FFFF))
                let sc = Int32(entry >> 32)
                if (col < low || col > high) && sc > score2 {
                    score2 = sc
                }
            }
        }

        return (bestScore, bestQEnd, bestTEnd, score2)
    }
}
