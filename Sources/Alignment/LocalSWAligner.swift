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
public struct LocalSWAligner: Sendable {

    /// Perform local Smith-Waterman alignment.
    ///
    /// - Parameters:
    ///   - query: Query sequence (2-bit encoded)
    ///   - target: Target/reference sequence (2-bit encoded)
    ///   - scoring: Scoring parameters
    /// - Returns: Alignment result with start/end positions, or nil if score <= 0
    public static func align(
        query: [UInt8], target: [UInt8], scoring: ScoringParameters
    ) -> LocalSWResult? {
        let qLen = query.count
        let tLen = target.count
        guard qLen > 0 && tLen > 0 else { return nil }

        let mat = scoring.scoringMatrix()
        let gapO = Int32(scoring.gapOpenPenalty)
        let gapE = Int32(scoring.gapExtendPenalty)

        // Forward pass: find best score and end positions
        var bestScore: Int32 = 0
        var bestQEnd: Int32 = -1
        var bestTEnd: Int32 = -1

        // DP arrays: H[j] = best score ending at target position j, E[j] = gap-in-query
        var H = [Int32](repeating: 0, count: tLen + 1)
        var E = [Int32](repeating: 0, count: tLen + 1)

        for i in 0..<qLen {
            let qBase = Int(min(query[i], 4))
            var f: Int32 = 0  // gap-in-target score
            var hPrev: Int32 = 0  // H[j-1] from previous column

            for j in 0..<tLen {
                let tBase = Int(min(target[j], 4))

                // Match/mismatch from diagonal
                let matchScore = mat[qBase * 5 + tBase]
                var h = hPrev + Int32(matchScore)

                // Gap in target (insertion in query): F
                f = max(H[j + 1] - gapO - gapE, f - gapE)

                // Gap in query (deletion in query): E
                E[j + 1] = max(H[j] - gapO - gapE, E[j + 1] - gapE)

                h = max(h, f)
                h = max(h, E[j + 1])
                h = max(h, 0)  // local alignment: no negative scores

                hPrev = H[j + 1]
                H[j + 1] = h

                if h > bestScore {
                    bestScore = h
                    bestQEnd = Int32(i)
                    bestTEnd = Int32(j)
                }
            }
            // Reset H[0] for next row (local alignment)
            hPrev = 0
        }

        guard bestScore > 0 else { return nil }

        // Reverse pass: find start positions by aligning reversed prefixes
        let revQuery = Array(query[0...Int(bestQEnd)].reversed())
        let revTarget = Array(target[0...Int(bestTEnd)].reversed())
        let revQLen = revQuery.count
        let revTLen = revTarget.count

        var bestRevScore: Int32 = 0
        var bestRevQEnd: Int32 = -1
        var bestRevTEnd: Int32 = -1

        var rH = [Int32](repeating: 0, count: revTLen + 1)
        var rE = [Int32](repeating: 0, count: revTLen + 1)

        for i in 0..<revQLen {
            let qBase = Int(min(revQuery[i], 4))
            var f: Int32 = 0
            var hPrev: Int32 = 0

            for j in 0..<revTLen {
                let tBase = Int(min(revTarget[j], 4))

                let matchScore = mat[qBase * 5 + tBase]
                var h = hPrev + Int32(matchScore)

                f = max(rH[j + 1] - gapO - gapE, f - gapE)
                rE[j + 1] = max(rH[j] - gapO - gapE, rE[j + 1] - gapE)

                h = max(h, f)
                h = max(h, rE[j + 1])
                h = max(h, 0)

                hPrev = rH[j + 1]
                rH[j + 1] = h

                if h > bestRevScore {
                    bestRevScore = h
                    bestRevQEnd = Int32(i)
                    bestRevTEnd = Int32(j)
                }
            }
        }

        // Convert reverse-pass "end" positions to original start positions
        let queryBegin = bestQEnd - bestRevQEnd
        let targetBegin = bestTEnd - bestRevTEnd

        return LocalSWResult(
            score: bestScore,
            queryBegin: queryBegin,
            queryEnd: bestQEnd,
            targetBegin: targetBegin,
            targetEnd: bestTEnd
        )
    }
}
