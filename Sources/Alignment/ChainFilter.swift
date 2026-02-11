import BWACore

/// Filters chains by weight, overlap, and marks primary/secondary.
/// Reimplements `mem_chain_flt()` from bwa-mem2's bwamem.cpp.
public struct ChainFilter: Sendable {

    /// Filter chains: remove low-weight chains and those significantly overlapping heavier chains.
    public static func filter(
        chains: inout [MemChain],
        scoring: ScoringParameters
    ) {
        guard !chains.isEmpty else { return }

        let dropRatio: Float = 0.50
        let minWeight = scoring.minSeedLength

        // Remove chains below minimum weight
        chains.removeAll { $0.weight < minWeight }
        guard !chains.isEmpty else { return }

        // Sort by weight descending
        chains.sort { $0.weight > $1.weight }

        // Mark kept/removed based on overlap with heavier chains
        for i in 0..<chains.count {
            if chains[i].kept == 0 { continue }

            let iBegin = chains[i].seeds.first?.qbeg ?? 0
            let iEnd = chains[i].seeds.last.map { $0.qbeg + $0.len } ?? 0

            for j in (i + 1)..<chains.count {
                if chains[j].kept == 0 { continue }

                let jBegin = chains[j].seeds.first?.qbeg ?? 0
                let jEnd = chains[j].seeds.last.map { $0.qbeg + $0.len } ?? 0

                // Calculate overlap
                let overlapBegin = max(iBegin, jBegin)
                let overlapEnd = min(iEnd, jEnd)

                if overlapBegin < overlapEnd {
                    let overlap = overlapEnd - overlapBegin
                    let jLen = jEnd - jBegin

                    if jLen > 0 && Float(overlap) / Float(jLen) > dropRatio &&
                       Float(chains[j].weight) < Float(chains[i].weight) * dropRatio {
                        chains[j].kept = 0
                    }
                }
            }
        }

        // Remove unkept chains
        chains.removeAll { $0.kept == 0 }
    }

    /// Mark secondary alignment regions based on overlap with primary regions.
    public static func markSecondary(
        regions: inout [MemAlnReg],
        maskLevel: Float
    ) {
        guard regions.count > 1 else { return }

        // Sort by score descending
        regions.sort { $0.score > $1.score }

        for i in 0..<regions.count {
            if regions[i].secondary >= 0 { continue }

            for j in (i + 1)..<regions.count {
                if regions[j].secondary >= 0 { continue }

                // Check overlap on query
                let qOverlapBegin = max(regions[i].qb, regions[j].qb)
                let qOverlapEnd = min(regions[i].qe, regions[j].qe)

                if qOverlapBegin < qOverlapEnd {
                    let overlap = Float(qOverlapEnd - qOverlapBegin)
                    let minLen = Float(min(regions[i].qe - regions[i].qb, regions[j].qe - regions[j].qb))

                    if minLen > 0 && overlap / minLen > maskLevel {
                        regions[j].secondary = Int32(i)
                    }
                }
            }
        }
    }
}
