import BWACore

/// Sorts, deduplicates, and patches (merges) overlapping alignment regions.
///
/// Reimplements `mem_sort_dedup_patch()` from bwa-mem2's `bwamem.cpp`.
/// Runs after chain extension but before secondary marking.
public struct RegionDedup: Sendable {

    private static let PATCH_MAX_R_BW: Float = 0.05
    private static let PATCH_MIN_SC_RATIO: Float = 0.90

    /// Sort, deduplicate, and patch overlapping alignment regions in-place.
    ///
    /// 1. Sort by (rid ASC, re ASC) to group same-chromosome regions
    /// 2. Remove redundant overlapping regions (>95% overlap on both ref and query)
    /// 3. Patch (merge) adjacent colinear regions via global re-alignment
    /// 4. Remove exact duplicates, sort by score descending
    ///
    /// - Parameters:
    ///   - regions: Alignment regions to process (modified in-place)
    ///   - query: 2-bit encoded read bases
    ///   - getReference: Closure to fetch reference subsequence at (position, length)
    ///   - genomeLength: Forward genome length (half of BWT total length)
    ///   - scoring: Scoring parameters
    public static func sortDedupPatch(
        regions: inout [MemAlnReg],
        query: [UInt8],
        getReference: (Int64, Int) -> [UInt8],
        genomeLength: Int64,
        scoring: ScoringParameters
    ) {
        guard regions.count > 1 else { return }

        // Step 1: Sort by (rid ASC, re ASC)
        regions.sort { a, b in
            if a.rid != b.rid { return a.rid < b.rid }
            return a.re < b.re
        }

        // Track n_comp per region (number of regions merged into this one)
        var nComp = [Int](repeating: 1, count: regions.count)

        // Step 2: Scan for redundant and patchable pairs
        let n = regions.count
        for i in 0..<n {
            // Skip deleted regions
            guard regions[i].qe > regions[i].qb else { continue }

            // Walk backwards through nearby regions on the same chromosome
            var j = i - 1
            while j >= 0 {
                // Must be same chromosome and within maxChainGap
                guard regions[j].rid == regions[i].rid else { break }
                guard regions[j].rb < regions[i].re + Int64(scoring.maxChainGap) else {
                    j -= 1
                    continue
                }
                // Skip deleted regions
                guard regions[j].qe > regions[j].qb else {
                    j -= 1
                    continue
                }

                // Check redundancy: >maskLevelRedun overlap on both ref and query
                let refOverlapBegin = max(regions[i].rb, regions[j].rb)
                let refOverlapEnd = min(regions[i].re, regions[j].re)
                let queryOverlapBegin = max(regions[i].qb, regions[j].qb)
                let queryOverlapEnd = min(regions[i].qe, regions[j].qe)

                if refOverlapBegin < refOverlapEnd && queryOverlapBegin < queryOverlapEnd {
                    let refOverlap = refOverlapBegin < refOverlapEnd
                        ? Float(refOverlapEnd - refOverlapBegin) : 0
                    let queryOverlap = queryOverlapBegin < queryOverlapEnd
                        ? Float(queryOverlapEnd - queryOverlapBegin) : 0

                    let minRefLen = Float(min(
                        regions[i].re - regions[i].rb,
                        regions[j].re - regions[j].rb
                    ))
                    let minQueryLen = Float(min(
                        regions[i].qe - regions[i].qb,
                        regions[j].qe - regions[j].qb
                    ))

                    if minRefLen > 0 && minQueryLen > 0
                        && refOverlap > scoring.maskLevelRedun * minRefLen
                        && queryOverlap > scoring.maskLevelRedun * minQueryLen
                    {
                        // Redundant: delete the lower-scoring one
                        if regions[i].score < regions[j].score {
                            regions[i].qe = regions[i].qb  // mark deleted
                            break  // i is deleted, stop scanning
                        } else {
                            regions[j].qe = regions[j].qb  // mark deleted
                            j -= 1
                            continue
                        }
                    }
                }

                // Check patchability: colinear regions that can be merged
                // q must be "left" of p in query coordinates (q.rb < p.rb)
                let (pIdx, qIdx): (Int, Int)
                if regions[j].rb < regions[i].rb {
                    pIdx = i; qIdx = j
                } else {
                    pIdx = j; qIdx = i
                }

                if regions[qIdx].qb < regions[pIdx].qb {
                    let mergedScore = patchRegion(
                        a: regions[qIdx],
                        b: regions[pIdx],
                        query: query,
                        getReference: getReference,
                        genomeLength: genomeLength,
                        scoring: scoring,
                        w: max(regions[qIdx].w, regions[pIdx].w)
                    )

                    if mergedScore > 0 {
                        // Merge: extend p to cover both, update score
                        regions[pIdx].rb = min(regions[pIdx].rb, regions[qIdx].rb)
                        regions[pIdx].re = max(regions[pIdx].re, regions[qIdx].re)
                        regions[pIdx].qb = min(regions[pIdx].qb, regions[qIdx].qb)
                        regions[pIdx].qe = max(regions[pIdx].qe, regions[qIdx].qe)
                        regions[pIdx].score = mergedScore
                        regions[pIdx].trueScore = mergedScore
                        regions[pIdx].w = max(regions[pIdx].w, regions[qIdx].w)
                        nComp[pIdx] += nComp[qIdx]
                        // Delete q
                        regions[qIdx].qe = regions[qIdx].qb
                        if qIdx == i { break }  // i was deleted
                    }
                }

                j -= 1
            }
        }

        // Step 3: Remove deleted entries, sort by score DESC
        regions.removeAll { $0.qe == $0.qb }
        regions.sort { $0.score > $1.score }

        // Step 4: Remove exact duplicates (same score, rb, qb)
        if regions.count > 1 {
            var i = 1
            while i < regions.count {
                if regions[i].score == regions[i - 1].score
                    && regions[i].rb == regions[i - 1].rb
                    && regions[i].qb == regions[i - 1].qb
                {
                    regions[i].qe = regions[i].qb  // mark deleted
                }
                i += 1
            }
            regions.removeAll { $0.qe == $0.qb }
        }
    }

    /// Attempt to patch (merge) two colinear regions via global re-alignment.
    ///
    /// - Parameters:
    ///   - a: Left region (lower query begin)
    ///   - b: Right region (higher query begin)
    ///   - query: 2-bit encoded read bases
    ///   - getReference: Closure to fetch reference at (position, length)
    ///   - genomeLength: Forward genome length
    ///   - scoring: Scoring parameters
    ///   - w: Maximum bandwidth from the two regions
    /// - Returns: Merged alignment score, or 0 if patch rejected
    private static func patchRegion(
        a: MemAlnReg,
        b: MemAlnReg,
        query: [UInt8],
        getReference: (Int64, Int) -> [UInt8],
        genomeLength: Int64,
        scoring: ScoringParameters,
        w: Int32
    ) -> Int32 {
        // Reject if different strands
        let aIsReverse = a.rb >= genomeLength
        let bIsReverse = b.rb >= genomeLength
        if aIsReverse != bIsReverse { return 0 }

        // Reject if not colinear: a must be before b in both query and ref
        guard a.qb < b.qb && a.qe < b.qe && a.re < b.re else { return 0 }

        // Compute bandwidth for the gap/overlap region
        let refGap = b.rb - a.re  // can be negative (overlap)
        let queryGap = Int64(b.qb) - Int64(a.qe)
        let bw = abs(Int32(truncatingIfNeeded: refGap - queryGap))
        let relBw: Float
        if a.re > b.rb {
            // Overlapping on reference
            relBw = Float(bw) / Float(a.re - b.rb)
        } else {
            relBw = Float(bw) / Float(b.rb - a.re)
        }

        // Apply thresholds: stricter for non-overlapping, more permissive for overlapping
        if a.re <= b.rb {
            // Non-overlapping (gap between regions)
            if bw > w << 1 || relBw >= PATCH_MAX_R_BW { return 0 }
        } else {
            // Overlapping on reference
            if bw > w << 2 || relBw >= PATCH_MAX_R_BW * 2 { return 0 }
        }

        // Perform global alignment over the merged span
        let mergedQb = a.qb
        let mergedQe = b.qe
        let mergedRb = a.rb
        let mergedRe = b.re

        guard mergedQe > mergedQb && mergedRe > mergedRb else { return 0 }

        let querySlice = Array(query[Int(mergedQb)..<Int(mergedQe)])
        let refSlice = getReference(mergedRb, Int(mergedRe - mergedRb))

        guard !querySlice.isEmpty && !refSlice.isEmpty else { return 0 }

        let mergedW = max(w, bw)
        let score: Int32 = querySlice.withUnsafeBufferPointer { qBuf in
            refSlice.withUnsafeBufferPointer { tBuf in
                GlobalAligner.align(
                    query: qBuf,
                    target: tBuf,
                    scoring: scoring,
                    w: mergedW
                ).score
            }
        }

        // Compute expected scores from query and reference lengths
        let queryLen = mergedQe - mergedQb
        let refLen = Int32(mergedRe - mergedRb)
        let expectedFromQuery = scoring.matchScore * queryLen
            - scoring.gapOpenPenalty
            - scoring.gapExtendPenalty * abs(queryLen - refLen)
        let expectedFromRef = scoring.matchScore * refLen
            - scoring.gapOpenPenalty
            - scoring.gapExtendPenalty * abs(queryLen - refLen)
        let expectedScore = max(expectedFromQuery, expectedFromRef)

        // Accept if score meets the minimum ratio threshold
        guard expectedScore > 0 else { return 0 }
        if Float(score) / Float(expectedScore) >= PATCH_MIN_SC_RATIO {
            return score
        }
        return 0
    }
}
