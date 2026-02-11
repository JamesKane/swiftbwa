import BWACore

/// Resolves paired-end alignments: mate rescue and proper pairing.
/// Reimplements paired-end logic from bwa-mem2's bwamem_pair.cpp.
public struct PairedEndResolver: Sendable {

    /// Pair two sets of alignment regions and resolve proper pairs.
    public static func resolve(
        regions1: [MemAlnReg],
        regions2: [MemAlnReg],
        insertStats: InsertSizeEstimator,
        scoring: ScoringParameters
    ) -> [(idx1: Int, idx2: Int, score: Int32)] {
        var pairs: [(idx1: Int, idx2: Int, score: Int32)] = []

        guard !regions1.isEmpty && !regions2.isEmpty else {
            return pairs
        }

        // Try all combinations of read1 x read2 alignments
        for (i, r1) in regions1.enumerated() {
            guard r1.secondary < 0 else { continue }  // Only primary

            for (j, r2) in regions2.enumerated() {
                guard r2.secondary < 0 else { continue }

                guard r1.rid == r2.rid else { continue }

                // Compute insert size
                let insertSize = abs(r2.re - r1.rb)

                // Score the pairing
                let pairScore = r1.score + r2.score
                let isProper = insertStats.isProperPair(Int64(insertSize))
                let bonus: Int32 = isProper ? scoring.unpairedPenalty : 0

                pairs.append((i, j, pairScore + bonus))
            }
        }

        // Sort by score descending
        pairs.sort { $0.score > $1.score }

        return pairs
    }
}
