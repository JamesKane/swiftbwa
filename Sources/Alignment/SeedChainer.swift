import BWACore

/// Converts SMEMs to positioned seeds and groups collinear seeds into chains.
/// Reimplements `mem_chain_seeds()` from bwa-mem2's bwamem.cpp.
public struct SeedChainer: Sendable {

    /// Build chains from SMEMs using suffix array lookups.
    ///
    /// - Parameters:
    ///   - smems: SMEMs found for one read
    ///   - getSAEntry: Closure that resolves SA interval position to reference coordinate
    ///   - metadata: Reference metadata for determining reference sequence IDs
    ///   - scoring: Scoring parameters
    ///   - readLength: Length of the read
    /// - Returns: Array of chains
    public static func chain(
        smems: [SMEM],
        getSAEntry: (Int64) -> Int64,
        metadata: ReferenceMetadata,
        scoring: ScoringParameters,
        readLength: Int32
    ) -> [MemChain] {
        let maxOcc = scoring.maxOccurrences
        let maxChainGap = scoring.maxChainGap

        var chains: [MemChain] = []

        for smem in smems {
            let seedLen = smem.queryEnd - smem.queryBegin
            guard seedLen >= scoring.minSeedLength else { continue }

            let occurrences = smem.count
            let step = occurrences > Int64(maxOcc) ? occurrences / Int64(maxOcc) : 1

            var saIdx: Int64 = 0
            var count: Int64 = 0
            while saIdx < occurrences && count < Int64(maxOcc) {
                let saPos = smem.k + saIdx
                let refPos = getSAEntry(saPos)

                // Determine reference sequence ID
                let genomeLen = metadata.totalLength
                var rbeg = refPos
                // Handle reverse complement: positions >= genomeLen are on reverse strand
                let isReverse = rbeg >= genomeLen
                if isReverse {
                    // Convert to forward coordinate for chaining
                    rbeg = genomeLen * 2 - 1 - rbeg - Int64(seedLen) + 1
                }

                let rid = metadata.sequenceID(for: min(rbeg, genomeLen - 1))
                let isAltContig = metadata.annotations[Int(rid)].isAlt

                let seed = MemSeed(
                    rbeg: refPos,
                    qbeg: smem.queryBegin,
                    len: seedLen,
                    score: seedLen * scoring.matchScore
                )

                // Try to merge into existing chain
                var merged = false
                for ci in stride(from: chains.count - 1, through: 0, by: -1) {
                    if chains[ci].rid != rid { continue }
                    if let lastSeed = chains[ci].seeds.last {
                        let refGap = abs(refPos - lastSeed.rbeg - Int64(lastSeed.len))
                        let queryGap = abs(Int64(seed.qbeg) - Int64(lastSeed.qbeg) - Int64(lastSeed.len))
                        let gap = max(refGap, queryGap)

                        if gap < Int64(maxChainGap) {
                            chains[ci].seeds.append(seed)
                            chains[ci].weight += seedLen
                            merged = true
                            break
                        }
                    }
                }

                if !merged {
                    var chain = MemChain()
                    chain.seeds = [seed]
                    chain.weight = seedLen
                    chain.rid = rid
                    chain.pos = refPos
                    chain.kept = 3
                    chain.isAlt = isAltContig
                    chains.append(chain)
                }

                saIdx += step
                count += 1
            }
        }

        // Sort seeds within each chain by reference position
        for i in 0..<chains.count {
            chains[i].seeds.sort { $0.rbeg < $1.rbeg }
            if let first = chains[i].seeds.first {
                chains[i].pos = first.rbeg
            }
        }

        return chains
    }
}
