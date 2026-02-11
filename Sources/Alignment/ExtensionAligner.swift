import BWACore

/// Drives Smith-Waterman extension left and right from seed chains.
/// Reimplements the extension logic from bwa-mem2's `mem_chain2aln_across_reads_V2()`.
public struct ExtensionAligner: Sendable {

    /// Extend a chain into an alignment region using banded Smith-Waterman.
    ///
    /// - Parameters:
    ///   - chain: The seed chain to extend
    ///   - read: The read sequence (2-bit encoded)
    ///   - getReference: Closure to get reference subsequence at (position, length)
    ///   - scoring: Scoring parameters
    /// - Returns: Array of alignment regions (typically one per chain, but may split)
    public static func extend(
        chain: MemChain,
        read: ReadSequence,
        getReference: (Int64, Int) -> [UInt8],
        scoring: ScoringParameters
    ) -> [MemAlnReg] {
        var regions: [MemAlnReg] = []
        let seeds = chain.seeds
        guard !seeds.isEmpty else { return regions }

        let readLen = Int32(read.length)
        let bandWidth = scoring.bandWidth

        // Process each seed as potential extension anchor
        var coveredRegions: [(qb: Int32, qe: Int32, rb: Int64, re: Int64)] = []

        for seedIdx in 0..<seeds.count {
            let seed = seeds[seedIdx]

            // Skip if this seed region is already covered by a previous extension
            let seedQEnd = seed.qbeg + seed.len
            let seedREnd = seed.rbeg + Int64(seed.len)

            var alreadyCovered = false
            for covered in coveredRegions {
                if seed.qbeg >= covered.qb && seedQEnd <= covered.qe {
                    alreadyCovered = true
                    break
                }
            }
            if alreadyCovered { continue }

            var reg = MemAlnReg()
            reg.w = bandWidth
            reg.rid = chain.rid
            reg.seedCov = seed.len
            reg.seedLen0 = Int32(seeds.count)

            // --- Left extension ---
            var leftScore: Int32 = 0
            var leftQLen: Int32 = 0
            var leftTLen: Int32 = 0

            if seed.qbeg > 0 {
                let queryLeft = Array(read.bases[0..<Int(seed.qbeg)].reversed())
                let targetLeftLen = min(Int(seed.qbeg) + Int(bandWidth), Int(seed.rbeg))

                if targetLeftLen > 0 {
                    let targetStart = seed.rbeg - Int64(targetLeftLen)
                    let targetLeft = getReference(targetStart, targetLeftLen).reversed()
                    let targetLeftArr = Array(targetLeft)

                    let result = queryLeft.withUnsafeBufferPointer { qBuf in
                        targetLeftArr.withUnsafeBufferPointer { tBuf in
                            BandedSWScalar.align(
                                query: qBuf,
                                target: tBuf,
                                scoring: scoring,
                                w: bandWidth,
                                h0: seed.score
                            )
                        }
                    }

                    leftScore = result.score
                    leftQLen = result.queryEnd
                    leftTLen = result.targetEnd
                }
            }

            // --- Right extension ---
            var rightScore: Int32 = 0
            var rightQLen: Int32 = 0
            var rightTLen: Int32 = 0

            if seedQEnd < readLen {
                let queryRight = Array(read.bases[Int(seedQEnd)..<read.length])
                let targetRightLen = min(Int(readLen - seedQEnd) + Int(bandWidth), 10000)

                if targetRightLen > 0 {
                    let targetRight = getReference(seedREnd, targetRightLen)

                    let result = queryRight.withUnsafeBufferPointer { qBuf in
                        targetRight.withUnsafeBufferPointer { tBuf in
                            BandedSWScalar.align(
                                query: qBuf,
                                target: tBuf,
                                scoring: scoring,
                                w: bandWidth,
                                h0: seed.score
                            )
                        }
                    }

                    rightScore = result.score
                    rightQLen = result.queryEnd
                    rightTLen = result.targetEnd
                }
            }

            // Compute final alignment coordinates
            reg.qb = seed.qbeg - leftQLen
            reg.qe = seedQEnd + rightQLen
            reg.rb = seed.rbeg - Int64(leftTLen)
            reg.re = seedREnd + Int64(rightTLen)
            reg.score = leftScore + Int32(seed.len) * scoring.matchScore + rightScore
            reg.trueScore = reg.score
            reg.sub = 0

            coveredRegions.append((reg.qb, reg.qe, reg.rb, reg.re))
            regions.append(reg)
        }

        return regions
    }
}
