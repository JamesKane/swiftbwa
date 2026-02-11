import BWACore

/// Generates CIGAR strings for alignment regions by running global DP on the
/// subsequences identified by Smith-Waterman extension.
///
/// Reimplements `bwa_gen_cigar2()` from bwa-mem2's `bwa.cpp` and the clipping
/// logic from `mem_reg2aln()` in `bwamem.cpp`.
public struct CIGARGenerator: Sendable {

    /// Generate a CIGAR string for an alignment region.
    ///
    /// The caller is responsible for extracting the correct reference subsequence
    /// and handling reverse-strand coordinate conversion. This function:
    /// 1. Runs GlobalAligner with bandwidth retry
    /// 2. Squeezes leading/trailing deletions (adjusting position)
    /// 3. Adds soft-clips based on query begin/end
    /// 4. Computes NM (edit distance)
    ///
    /// - Parameters:
    ///   - querySegment: 2-bit encoded query bases for [qb..qe) (already RC'd if reverse)
    ///   - refSegment: 2-bit encoded reference bases for [rb..re) (forward strand)
    ///   - qb: Query begin in read coordinates
    ///   - qe: Query end in read coordinates
    ///   - readLength: Total read length (for soft-clip computation)
    ///   - isReverse: Whether alignment is on reverse strand
    ///   - trueScore: Expected alignment score (for retry loop)
    ///   - initialW: Initial band width
    ///   - scoring: Scoring parameters
    ///   - refPos: Reference start position (may be adjusted by leading deletion squeeze)
    public static func generate(
        querySegment: [UInt8],
        refSegment: [UInt8],
        qb: Int32,
        qe: Int32,
        readLength: Int,
        isReverse: Bool,
        trueScore: Int32,
        initialW: Int32,
        scoring: ScoringParameters,
        refPos: Int64
    ) -> CIGARResult {
        // Infer initial band width from score and lengths
        var w = inferBandWidth(
            queryLen: qe - qb,
            refLen: Int32(refSegment.count),
            score: trueScore,
            scoring: scoring
        )
        w = max(w, initialW)

        // Run global alignment with retry loop (matching bwa-mem2's retry at bwamem.cpp:1758-1766)
        var result: GlobalAlignResult
        var retries = 0
        let maxRetries = 3

        result = querySegment.withUnsafeBufferPointer { qBuf in
            refSegment.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: w)
            }
        }

        while result.score < trueScore && retries < maxRetries {
            w *= 2
            retries += 1
            result = querySegment.withUnsafeBufferPointer { qBuf in
                refSegment.withUnsafeBufferPointer { tBuf in
                    GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: w)
                }
            }
        }

        var cigar = result.cigar
        var pos = refPos

        // Squeeze leading deletions (adjust position instead)
        while !cigar.isEmpty {
            let op = cigar[0] & 0xF
            if op == CIGAROp.deletion.rawValue {
                let len = Int64(cigar[0] >> 4)
                pos += len
                cigar.removeFirst()
            } else {
                break
            }
        }

        // Squeeze trailing deletions
        while !cigar.isEmpty {
            let op = cigar[cigar.count - 1] & 0xF
            if op == CIGAROp.deletion.rawValue {
                cigar.removeLast()
            } else {
                break
            }
        }

        // Compute NM from the (possibly trimmed) CIGAR
        // NM is computed on the aligned portion after deletion squeeze
        let nm = computeNMFromCigar(cigar: cigar, query: querySegment, target: refSegment)

        // Add soft-clips
        let intQb = Int(qb)
        let intQe = Int(qe)
        if isReverse {
            let clip5 = readLength - intQe
            let clip3 = intQb
            if clip5 > 0 {
                cigar.insert(UInt32(clip5) << 4 | CIGAROp.softClip.rawValue, at: 0)
            }
            if clip3 > 0 {
                cigar.append(UInt32(clip3) << 4 | CIGAROp.softClip.rawValue)
            }
        } else {
            if intQb > 0 {
                cigar.insert(UInt32(intQb) << 4 | CIGAROp.softClip.rawValue, at: 0)
            }
            if intQe < readLength {
                cigar.append(UInt32(readLength - intQe) << 4 | CIGAROp.softClip.rawValue)
            }
        }

        return CIGARResult(cigar: cigar, nm: nm, score: result.score, pos: pos)
    }

    /// Infer band width from score difference and gap penalties.
    ///
    /// Matches the logic in bwa-mem2's `mem_reg2aln()` (bwamem.cpp:1742-1755).
    static func inferBandWidth(
        queryLen: Int32,
        refLen: Int32,
        score: Int32,
        scoring: ScoringParameters
    ) -> Int32 {
        let maxLen = max(queryLen, refLen)
        let minLen = min(queryLen, refLen)
        let diff = maxLen - minLen

        let expectedScore = minLen * scoring.matchScore

        if score >= expectedScore {
            return diff + 3
        }

        let deficit = expectedScore - score
        let gapCost = scoring.gapOpenPenalty + scoring.gapExtendPenalty
        let matchMismatchDelta = scoring.matchScore + scoring.mismatchPenalty

        let estimatedErrors = deficit / min(gapCost, matchMismatchDelta)

        return max(diff + 3, estimatedErrors)
    }

    /// Compute NM from CIGAR by comparing query and target bases.
    private static func computeNMFromCigar(
        cigar: [UInt32],
        query: [UInt8],
        target: [UInt8]
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
            default:
                break
            }
        }

        return nm
    }
}
