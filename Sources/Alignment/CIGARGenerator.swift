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
        refPos: Int64,
        scoringMatrix: [Int8]? = nil
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

        let mat = scoringMatrix ?? scoring.scoringMatrix()

        // Fast path: when query and ref are the same length and the score deficit
        // is too small for any indel pair (open+extend for both ins and del), the
        // alignment must be pure matches/mismatches — CIGAR is just "{len}M".
        // This skips the entire GlobalAligner DP for ~90% of Illumina reads.
        let queryLen = querySegment.count
        let refLen = refSegment.count
        let minIndelCost = scoring.gapOpenPenalty + scoring.gapExtendPenalty
            + scoring.gapOpenPenaltyDeletion + scoring.gapExtendPenaltyDeletion
        let deficit = Int32(queryLen) * scoring.matchScore - trueScore

        if queryLen == refLen && queryLen > 0 && deficit < minIndelCost {
            let cigarOp = UInt32(queryLen) << 4 | CIGAROp.match.rawValue
            var cigar = [cigarOp]
            let nm = computeNMFromCigar(cigar: cigar, query: querySegment,
                                         target: refSegment)
            let md = generateMD(cigar: cigar, query: querySegment,
                                target: refSegment)

            addSoftClips(&cigar, qb: qb, qe: qe, readLength: readLength, isReverse: isReverse)

            return CIGARResult(cigar: cigar, nm: nm, md: md, score: trueScore, pos: refPos)
        }

        result = querySegment.withUnsafeBufferPointer { qBuf in
            refSegment.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: w, scoringMatrix: mat)
            }
        }

        while result.score < trueScore && retries < maxRetries {
            w *= 2
            retries += 1
            result = querySegment.withUnsafeBufferPointer { qBuf in
                refSegment.withUnsafeBufferPointer { tBuf in
                    GlobalAligner.align(query: qBuf, target: tBuf, scoring: scoring, w: w, scoringMatrix: mat)
                }
            }
        }

        var cigar = result.cigar
        var pos = refPos
        let (posAdj, refOffset) = squeezeDeletions(&cigar)
        pos += posAdj

        let nm = computeNMFromCigar(cigar: cigar, query: querySegment,
                                     target: refSegment, targetOffset: refOffset)
        let md = generateMD(cigar: cigar, query: querySegment,
                            target: refSegment, targetOffset: refOffset)

        addSoftClips(&cigar, qb: qb, qe: qe, readLength: readLength, isReverse: isReverse)

        return CIGARResult(cigar: cigar, nm: nm, md: md, score: result.score, pos: pos)
    }

    /// Post-process a GPU-produced CIGAR: squeeze leading/trailing deletions,
    /// compute NM/MD, add soft-clips. Used by the GPU global SW path.
    ///
    /// The raw CIGAR comes from the GPU `global_sw` kernel (banded NW + traceback).
    /// This applies the same post-processing as `generate()` lines 117-172.
    public static func postProcess(
        cigar rawCigar: [UInt32],
        score: Int32,
        querySegment: [UInt8],
        refSegment: [UInt8],
        qb: Int32,
        qe: Int32,
        readLength: Int,
        isReverse: Bool,
        refPos: Int64
    ) -> CIGARResult {
        var cigar = rawCigar
        var pos = refPos
        let (posAdj, refOffset) = squeezeDeletions(&cigar)
        pos += posAdj

        let nm = computeNMFromCigar(cigar: cigar, query: querySegment,
                                     target: refSegment, targetOffset: refOffset)
        let md = generateMD(cigar: cigar, query: querySegment,
                            target: refSegment, targetOffset: refOffset)

        addSoftClips(&cigar, qb: qb, qe: qe, readLength: readLength, isReverse: isReverse)

        return CIGARResult(cigar: cigar, nm: nm, md: md, score: score, pos: pos)
    }

    /// Infer band width from score difference and gap penalties.
    ///
    /// Matches the logic in bwa-mem2's `mem_reg2aln()` (bwamem.cpp:1742-1755).
    public static func inferBandWidth(
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
        target: [UInt8],
        targetOffset: Int = 0
    ) -> Int32 {
        var nm: Int32 = 0
        var qi = 0
        var ti = targetOffset

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

    /// Generate the MD tag string from CIGAR and reference bases.
    ///
    /// MD format: match counts interspersed with mismatched ref bases and
    /// `^`-prefixed deleted ref bases. E.g. `50A^GC49` means 50 matches,
    /// mismatch (ref=A), 2bp deletion (ref=GC), 49 matches.
    ///
    /// Only M/I/D operations affect MD. Soft-clips are skipped.
    private static func generateMD(
        cigar: [UInt32],
        query: [UInt8],
        target: [UInt8],
        targetOffset: Int = 0
    ) -> String {
        var md = ""
        var matchCount = 0
        var qi = 0
        var ti = targetOffset

        for c in cigar {
            let op = c & 0xF
            let len = Int(c >> 4)

            switch op {
            case CIGAROp.match.rawValue:
                for _ in 0..<len {
                    if qi < query.count && ti < target.count {
                        if query[qi] == target[ti] {
                            matchCount += 1
                        } else {
                            md += String(matchCount)
                            md += baseChar(target[ti])
                            matchCount = 0
                        }
                        qi += 1
                        ti += 1
                    }
                }
            case CIGAROp.insertion.rawValue:
                qi += len  // insertions don't appear in MD
            case CIGAROp.deletion.rawValue:
                md += String(matchCount)
                md += "^"
                for _ in 0..<len {
                    if ti < target.count {
                        md += baseChar(target[ti])
                        ti += 1
                    }
                }
                matchCount = 0
            default:
                break
            }
        }

        md += String(matchCount)
        return md
    }

    /// Remove leading and trailing deletion ops from CIGAR, adjusting position.
    /// Returns (position adjustment from leading deletions, ref bases consumed by leading deletions).
    private static func squeezeDeletions(_ cigar: inout [UInt32]) -> (posAdj: Int64, refOffset: Int) {
        var posAdj: Int64 = 0
        var refOffset = 0
        while !cigar.isEmpty {
            let op = cigar[0] & 0xF
            if op == CIGAROp.deletion.rawValue {
                let len = Int(cigar[0] >> 4)
                posAdj += Int64(len)
                refOffset += len
                cigar.removeFirst()
            } else {
                break
            }
        }
        while !cigar.isEmpty {
            let op = cigar[cigar.count - 1] & 0xF
            if op == CIGAROp.deletion.rawValue {
                cigar.removeLast()
            } else {
                break
            }
        }
        return (posAdj, refOffset)
    }

    /// Add soft-clip ops to the front/back of CIGAR based on alignment bounds.
    private static func addSoftClips(
        _ cigar: inout [UInt32], qb: Int32, qe: Int32, readLength: Int, isReverse: Bool
    ) {
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
    }

    /// Convert 2-bit encoded base to character.
    private static func baseChar(_ base: UInt8) -> String {
        switch base {
        case 0: return "A"
        case 1: return "C"
        case 2: return "G"
        case 3: return "T"
        default: return "N"
        }
    }
}
