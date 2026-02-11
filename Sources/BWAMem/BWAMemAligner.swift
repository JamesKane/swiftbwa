import BWACore
import FMIndex
import Alignment
import Htslib

/// Main BWA-MEM alignment pipeline.
/// Orchestrates SMEM finding, chaining, extension, and SAM output.
public actor BWAMemAligner {
    public let index: FMIndex
    public let options: BWAMemOptions

    public init(index: FMIndex, options: BWAMemOptions = BWAMemOptions()) {
        self.index = index
        self.options = options
    }

    /// Align a single read against the reference.
    public func alignRead(_ read: ReadSequence) -> [MemAlnReg] {
        let scoring = options.scoring

        // Phase 1: Find SMEMs
        let smems = SMEMFinder.findAllSMEMs(
            query: read.bases,
            bwt: index.bwt,
            minSeedLen: scoring.minSeedLength
        )

        guard !smems.isEmpty else { return [] }

        // Phase 2: Chain seeds
        var chains = SeedChainer.chain(
            smems: smems,
            getSAEntry: { [index] pos in
                index.suffixArray.entry(at: pos)
            },
            metadata: index.metadata,
            scoring: scoring,
            readLength: Int32(read.length)
        )

        // Phase 3: Filter chains
        ChainFilter.filter(chains: &chains, scoring: scoring)

        guard !chains.isEmpty else { return [] }

        // Phase 4: Extend chains with Smith-Waterman
        var regions: [MemAlnReg] = []

        for chain in chains {
            let chainRegions = ExtensionAligner.extend(
                chain: chain,
                read: read,
                getReference: { [index] pos, length in
                    let safeLen = min(length, Int(index.packedRef.length - pos))
                    guard safeLen > 0 && pos >= 0 else { return [] }
                    return index.packedRef.subsequence(from: pos, length: safeLen)
                },
                scoring: scoring
            )
            regions.append(contentsOf: chainRegions)
        }

        // Phase 5: Mark secondary alignments
        ChainFilter.markSecondary(regions: &regions, maskLevel: scoring.maskLevel)

        return regions
    }

    /// Align a batch of reads and write SAM/BAM output.
    ///
    /// Reads are aligned in parallel using a task group, then output is written
    /// sequentially in input order. The `outputFile` parameter is borrowed because
    /// `HTSFile` is a move-only type.
    public func alignBatch(
        reads: [ReadSequence],
        outputFile: borrowing HTSFile,
        header: SAMHeader
    ) async throws {
        // Process reads in parallel using TaskGroup
        let results = await withTaskGroup(
            of: (Int, [MemAlnReg]).self,
            returning: [(Int, [MemAlnReg])].self
        ) { group in
            for (idx, read) in reads.enumerated() {
                group.addTask { [self] in
                    let regions = await self.alignRead(read)
                    return (idx, regions)
                }
            }

            var collected: [(Int, [MemAlnReg])] = []
            for await result in group {
                collected.append(result)
            }
            return collected.sorted { $0.0 < $1.0 }
        }

        // Write output sequentially (ordered by input)
        for (idx, regions) in results {
            let read = reads[idx]

            if regions.isEmpty {
                let record = try SAMOutputBuilder.buildUnmappedRecord(read: read)
                try outputFile.write(record: record, header: header)
            } else {
                for (regIdx, region) in regions.enumerated() {
                    guard region.score >= options.scoring.minOutputScore else { continue }
                    let isPrimary = region.secondary < 0 && regIdx == 0

                    let genomeLen = index.genomeLength
                    let isReverse = region.rb >= genomeLen

                    // Extract query and reference segments for global alignment
                    let qb = Int(region.qb)
                    let qe = Int(region.qe)

                    let querySegment: [UInt8]
                    if isReverse {
                        let rc = read.reverseComplement()
                        let rcQb = read.length - qe
                        let rcQe = read.length - qb
                        querySegment = Array(rc[rcQb..<rcQe])
                    } else {
                        querySegment = Array(read.bases[qb..<qe])
                    }

                    var rb = region.rb
                    var re = region.re
                    if isReverse {
                        let fwdRb = 2 * genomeLen - re
                        let fwdRe = 2 * genomeLen - rb
                        rb = fwdRb
                        re = fwdRe
                    }
                    let refLen = Int(re - rb)
                    let safeRefLen = min(refLen, Int(index.packedRef.length - rb))
                    let refSegment: [UInt8]
                    if safeRefLen > 0 && rb >= 0 {
                        refSegment = index.packedRef.subsequence(from: rb, length: safeRefLen)
                    } else {
                        refSegment = []
                    }

                    let cigarResult = CIGARGenerator.generate(
                        querySegment: querySegment,
                        refSegment: refSegment,
                        qb: region.qb,
                        qe: region.qe,
                        readLength: read.length,
                        isReverse: isReverse,
                        trueScore: region.trueScore,
                        initialW: region.w,
                        scoring: options.scoring,
                        refPos: rb
                    )
                    let cigar = cigarResult.cigar
                    let nm = cigarResult.nm

                    let record = try SAMOutputBuilder.buildRecord(
                        read: read,
                        region: region,
                        allRegions: regions,
                        metadata: index.metadata,
                        scoring: options.scoring,
                        cigar: cigar,
                        nm: nm,
                        md: cigarResult.md,
                        isPrimary: isPrimary,
                        adjustedPos: cigarResult.pos
                    )
                    try outputFile.write(record: record, header: header)
                }
            }
        }
    }
}
