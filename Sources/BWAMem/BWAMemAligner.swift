import BWACore
import FMIndex
import Alignment
import Htslib
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif canImport(Musl)
import Musl
#endif

/// Result of CIGAR generation for a single alignment.
struct CIGARInfo: Sendable {
    var cigar: [UInt32]
    var nm: Int32
    var md: String
    var pos: Int64
    var isReverse: Bool
    /// Reference bases consumed by this alignment's CIGAR.
    var refConsumed: Int64
    /// CIGAR as a SAM-format string (e.g., "50M2I48M").
    var cigarString: String
}

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

    /// Generate CIGAR for a single alignment region.
    func generateCIGAR(read: ReadSequence, region: MemAlnReg) -> CIGARInfo {
        let genomeLen = index.genomeLength
        let isReverse = region.rb >= genomeLen

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

        // Compute reference consumed and CIGAR string
        var refConsumed: Int64 = 0
        var cigarStr = ""
        for packed in cigarResult.cigar {
            let op = CIGAROperation(rawValue: packed)
            let len = op.length
            if op.consumesReference {
                refConsumed += Int64(len)
            }
            cigarStr += "\(len)\(op.character)"
        }

        return CIGARInfo(
            cigar: cigarResult.cigar,
            nm: cigarResult.nm,
            md: cigarResult.md,
            pos: cigarResult.pos,
            isReverse: isReverse,
            refConsumed: refConsumed,
            cigarString: cigarStr
        )
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

                    let cigarInfo = generateCIGAR(read: read, region: region)

                    let record = try SAMOutputBuilder.buildRecord(
                        read: read,
                        region: region,
                        allRegions: regions,
                        metadata: index.metadata,
                        scoring: options.scoring,
                        cigar: cigarInfo.cigar,
                        nm: cigarInfo.nm,
                        md: cigarInfo.md,
                        isPrimary: isPrimary,
                        adjustedPos: cigarInfo.pos
                    )
                    try outputFile.write(record: record, header: header)
                }
            }
        }
    }

    /// Align paired-end reads and write SAM/BAM output.
    ///
    /// Both ends are aligned in parallel, insert size is estimated from the first
    /// batch, pairs are resolved using z-score-based scoring, and output is written
    /// with full paired-end flags.
    public func alignPairedBatch(
        reads1: [ReadSequence],
        reads2: [ReadSequence],
        outputFile: borrowing HTSFile,
        header: SAMHeader
    ) async throws {
        let pairCount = min(reads1.count, reads2.count)
        let genomeLen = index.genomeLength

        // Phase 1: Align all reads in parallel
        let allRegions1 = await alignAllReads(reads1)
        let allRegions2 = await alignAllReads(reads2)

        // Phase 2: Estimate insert size distribution
        let dist = InsertSizeEstimator.estimate(
            regions1: allRegions1,
            regions2: allRegions2,
            genomeLength: genomeLen
        )

        let primaryStats = dist.stats[dist.primaryOrientation.rawValue]
        if !primaryStats.failed {
            fputs("[PE] Insert size: mean=\(String(format: "%.1f", primaryStats.mean)), "
                  + "stddev=\(String(format: "%.1f", primaryStats.stddev)), "
                  + "orientation=\(dist.primaryOrientation), "
                  + "n=\(primaryStats.count)\n", stderr)
        } else {
            fputs("[PE] Warning: insert size estimation failed, "
                  + "treating as unpaired for scoring\n", stderr)
        }

        // Phase 3: For each pair, resolve and write output
        for i in 0..<pairCount {
            let read1 = reads1[i]
            let read2 = reads2[i]
            let regions1 = allRegions1[i]
            let regions2 = allRegions2[i]

            let r1Empty = regions1.isEmpty
            let r2Empty = regions2.isEmpty

            if r1Empty && r2Empty {
                // Both unmapped
                let pe1 = PairedEndInfo(
                    isRead1: true, isProperPair: false,
                    mateTid: -1, matePos: -1,
                    mateIsReverse: false, mateIsUnmapped: true,
                    tlen: 0
                )
                let pe2 = PairedEndInfo(
                    isRead1: false, isProperPair: false,
                    mateTid: -1, matePos: -1,
                    mateIsReverse: false, mateIsUnmapped: true,
                    tlen: 0
                )
                let rec1 = try SAMOutputBuilder.buildUnmappedRecord(read: read1, pairedEnd: pe1)
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildUnmappedRecord(read: read2, pairedEnd: pe2)
                try outputFile.write(record: rec2, header: header)
                continue
            }

            // Try to resolve best pair
            let decision = PairedEndResolver.resolve(
                regions1: regions1,
                regions2: regions2,
                dist: dist,
                genomeLength: genomeLen,
                scoring: options.scoring
            )

            if let decision = decision {
                // Both mapped and paired
                let reg1 = regions1[decision.idx1]
                let reg2 = regions2[decision.idx2]

                let cigar1 = generateCIGAR(read: read1, region: reg1)
                let cigar2 = generateCIGAR(read: read2, region: reg2)

                let (rid1, localPos1) = index.metadata.decodePosition(cigar1.pos)
                let (rid2, localPos2) = index.metadata.decodePosition(cigar2.pos)

                let (tlen1, tlen2) = PairedEndResolver.computeTLEN(
                    pos1: localPos1, isReverse1: cigar1.isReverse, refLen1: cigar1.refConsumed,
                    pos2: localPos2, isReverse2: cigar2.isReverse, refLen2: cigar2.refConsumed
                )

                let pe1 = PairedEndInfo(
                    isRead1: true, isProperPair: decision.isProperPair,
                    mateTid: rid2, matePos: localPos2,
                    mateIsReverse: cigar2.isReverse, mateIsUnmapped: false,
                    tlen: tlen1, mateCigarString: cigar2.cigarString
                )
                let pe2 = PairedEndInfo(
                    isRead1: false, isProperPair: decision.isProperPair,
                    mateTid: rid1, matePos: localPos1,
                    mateIsReverse: cigar1.isReverse, mateIsUnmapped: false,
                    tlen: tlen2, mateCigarString: cigar1.cigarString
                )

                let rec1 = try SAMOutputBuilder.buildRecord(
                    read: read1, region: reg1, allRegions: regions1,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar1.cigar, nm: cigar1.nm, md: cigar1.md,
                    isPrimary: true, adjustedPos: cigar1.pos, pairedEnd: pe1
                )
                try outputFile.write(record: rec1, header: header)

                let rec2 = try SAMOutputBuilder.buildRecord(
                    read: read2, region: reg2, allRegions: regions2,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar2.cigar, nm: cigar2.nm, md: cigar2.md,
                    isPrimary: true, adjustedPos: cigar2.pos, pairedEnd: pe2
                )
                try outputFile.write(record: rec2, header: header)
            } else if !r1Empty && r2Empty {
                // Read 1 mapped, read 2 unmapped
                writeMappedUnmappedPair(
                    mappedRead: read1, mappedRegions: regions1,
                    unmappedRead: read2,
                    mappedIsRead1: true,
                    outputFile: outputFile, header: header
                )
            } else if r1Empty && !r2Empty {
                // Read 1 unmapped, read 2 mapped
                writeMappedUnmappedPair(
                    mappedRead: read2, mappedRegions: regions2,
                    unmappedRead: read1,
                    mappedIsRead1: false,
                    outputFile: outputFile, header: header
                )
            } else {
                // Both have regions but no valid pairing (e.g., different chromosomes)
                // Write best region for each as unpaired
                let bestReg1 = regions1[0]
                let bestReg2 = regions2[0]
                let cigar1 = generateCIGAR(read: read1, region: bestReg1)
                let cigar2 = generateCIGAR(read: read2, region: bestReg2)

                let (rid1, localPos1) = index.metadata.decodePosition(cigar1.pos)
                let (rid2, localPos2) = index.metadata.decodePosition(cigar2.pos)

                let pe1 = PairedEndInfo(
                    isRead1: true, isProperPair: false,
                    mateTid: rid2, matePos: localPos2,
                    mateIsReverse: cigar2.isReverse, mateIsUnmapped: false,
                    tlen: 0, mateCigarString: cigar2.cigarString
                )
                let pe2 = PairedEndInfo(
                    isRead1: false, isProperPair: false,
                    mateTid: rid1, matePos: localPos1,
                    mateIsReverse: cigar1.isReverse, mateIsUnmapped: false,
                    tlen: 0, mateCigarString: cigar1.cigarString
                )

                let rec1 = try SAMOutputBuilder.buildRecord(
                    read: read1, region: bestReg1, allRegions: regions1,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar1.cigar, nm: cigar1.nm, md: cigar1.md,
                    isPrimary: true, adjustedPos: cigar1.pos, pairedEnd: pe1
                )
                try outputFile.write(record: rec1, header: header)

                let rec2 = try SAMOutputBuilder.buildRecord(
                    read: read2, region: bestReg2, allRegions: regions2,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar2.cigar, nm: cigar2.nm, md: cigar2.md,
                    isPrimary: true, adjustedPos: cigar2.pos, pairedEnd: pe2
                )
                try outputFile.write(record: rec2, header: header)
            }
        }
    }

    // MARK: - Private Helpers

    /// Align all reads in parallel and return regions indexed by read position.
    private func alignAllReads(_ reads: [ReadSequence]) async -> [[MemAlnReg]] {
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
        return results.map { $0.1 }
    }

    /// Write a pair where one read is mapped and the other is unmapped.
    private func writeMappedUnmappedPair(
        mappedRead: ReadSequence,
        mappedRegions: [MemAlnReg],
        unmappedRead: ReadSequence,
        mappedIsRead1: Bool,
        outputFile: borrowing HTSFile,
        header: SAMHeader
    ) {
        let region = mappedRegions[0]
        let cigar = generateCIGAR(read: mappedRead, region: region)
        let (rid, localPos) = index.metadata.decodePosition(cigar.pos)

        // Mapped read's PE info: mate is unmapped
        let mappedPE = PairedEndInfo(
            isRead1: mappedIsRead1, isProperPair: false,
            mateTid: rid, matePos: localPos,
            mateIsReverse: false, mateIsUnmapped: true,
            tlen: 0
        )

        // Unmapped read's PE info: mate is mapped
        let unmappedPE = PairedEndInfo(
            isRead1: !mappedIsRead1, isProperPair: false,
            mateTid: rid, matePos: localPos,
            mateIsReverse: cigar.isReverse, mateIsUnmapped: false,
            tlen: 0
        )

        do {
            if mappedIsRead1 {
                let rec1 = try SAMOutputBuilder.buildRecord(
                    read: mappedRead, region: region, allRegions: mappedRegions,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar.cigar, nm: cigar.nm, md: cigar.md,
                    isPrimary: true, adjustedPos: cigar.pos, pairedEnd: mappedPE
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: unmappedRead, pairedEnd: unmappedPE
                )
                try outputFile.write(record: rec2, header: header)
            } else {
                let rec1 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: unmappedRead, pairedEnd: unmappedPE
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildRecord(
                    read: mappedRead, region: region, allRegions: mappedRegions,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar.cigar, nm: cigar.nm, md: cigar.md,
                    isPrimary: true, adjustedPos: cigar.pos, pairedEnd: mappedPE
                )
                try outputFile.write(record: rec2, header: header)
            }
        } catch {
            fputs("[PE] Error writing mapped/unmapped pair: \(error)\n", stderr)
        }
    }
}
