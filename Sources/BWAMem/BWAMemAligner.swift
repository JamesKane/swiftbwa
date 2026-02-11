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

/// Classification of an alignment segment for output.
struct AlnSegment: Sendable {
    let regionIndex: Int
    let cigarInfo: CIGARInfo
    var mapq: UInt8
    let isPrimary: Bool
    let isSupplementary: Bool
    let isSecondary: Bool
    let rname: String
    let localPos: Int64
    let rid: Int32
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
    nonisolated public func alignRead(_ read: ReadSequence) -> [MemAlnReg] {
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

        // Safety net: ensure isAlt is set from metadata (matches bwa-mem2 line 1166-1167)
        for i in 0..<regions.count {
            if regions[i].rid >= 0
                && regions[i].rid < index.metadata.numSequences
                && index.metadata.annotations[Int(regions[i].rid)].isAlt {
                regions[i].isAlt = true
            }
        }

        // Phase 5: Mark secondary alignments (ALT-aware if any ALT regions)
        let hasAlt = regions.contains { $0.isAlt }
        if hasAlt {
            ChainFilter.markSecondaryALT(
                regions: &regions, maskLevel: scoring.maskLevel, scoring: scoring
            )
        } else {
            ChainFilter.markSecondary(regions: &regions, maskLevel: scoring.maskLevel)
        }

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
        // Process reads in parallel using TaskGroup with sliding-window concurrency
        let maxConcurrency = options.scoring.numThreads

        let results = await withTaskGroup(
            of: (Int, [MemAlnReg]).self,
            returning: [(Int, [MemAlnReg])].self
        ) { group in
            var nextIdx = 0
            var collected: [(Int, [MemAlnReg])] = []

            // Seed initial batch
            while nextIdx < min(maxConcurrency, reads.count) {
                let idx = nextIdx
                let read = reads[idx]
                group.addTask { [self] in
                    (idx, self.alignRead(read))
                }
                nextIdx += 1
            }

            // As each completes, launch next
            for await result in group {
                collected.append(result)
                if nextIdx < reads.count {
                    let idx = nextIdx
                    let read = reads[idx]
                    group.addTask { [self] in
                        (idx, self.alignRead(read))
                    }
                    nextIdx += 1
                }
            }

            return collected.sorted { $0.0 < $1.0 }
        }

        // Write output sequentially (ordered by input)
        let rgID = options.readGroupID
        for (idx, regions) in results {
            let read = reads[idx]

            if regions.isEmpty {
                let record = try SAMOutputBuilder.buildUnmappedRecord(
                    read: read, readGroupID: rgID
                )
                try outputFile.write(record: record, header: header)
            } else {
                try emitSingleEndAlignments(
                    read: read,
                    regions: regions,
                    outputFile: outputFile,
                    header: header
                )
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

        // Phase 2.5: Mate rescue
        var mutableRegions1 = allRegions1
        var mutableRegions2 = allRegions2

        if !primaryStats.failed {
            let scoring = options.scoring
            for i in 0..<pairCount {
                // Rescue read2 using read1's alignments as templates
                let templates1 = selectRescueCandidates(mutableRegions1[i], scoring: scoring)
                let rescued2 = MateRescue.rescue(
                    templateRegions: templates1,
                    mateRead: reads2[i],
                    mateRegions: mutableRegions2[i],
                    dist: dist,
                    genomeLength: genomeLen,
                    packedRef: index.packedRef,
                    metadata: index.metadata,
                    scoring: scoring
                )
                mutableRegions2[i].append(contentsOf: rescued2)

                // Rescue read1 using read2's alignments (symmetric)
                let templates2 = selectRescueCandidates(mutableRegions2[i], scoring: scoring)
                let rescued1 = MateRescue.rescue(
                    templateRegions: templates2,
                    mateRead: reads1[i],
                    mateRegions: mutableRegions1[i],
                    dist: dist,
                    genomeLength: genomeLen,
                    packedRef: index.packedRef,
                    metadata: index.metadata,
                    scoring: scoring
                )
                mutableRegions1[i].append(contentsOf: rescued1)

                // Re-mark secondaries after adding rescued regions
                if mutableRegions1[i].count > 1 {
                    ChainFilter.markSecondary(
                        regions: &mutableRegions1[i], maskLevel: scoring.maskLevel
                    )
                }
                if mutableRegions2[i].count > 1 {
                    ChainFilter.markSecondary(
                        regions: &mutableRegions2[i], maskLevel: scoring.maskLevel
                    )
                }
            }
        }

        // Phase 3: For each pair, resolve and write output
        let rgID = options.readGroupID
        for i in 0..<pairCount {
            let read1 = reads1[i]
            let read2 = reads2[i]
            let regions1 = mutableRegions1[i]
            let regions2 = mutableRegions2[i]

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
                let rec1 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: read1, pairedEnd: pe1, readGroupID: rgID
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: read2, pairedEnd: pe2, readGroupID: rgID
                )
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

                // Emit primary read1 using emitSingleEndAlignments for full supplementary handling
                try emitSingleEndAlignments(
                    read: read1, regions: regions1,
                    outputFile: outputFile, header: header,
                    pairedEnd: pe1
                )

                // Emit primary read2 using emitSingleEndAlignments for full supplementary handling
                try emitSingleEndAlignments(
                    read: read2, regions: regions2,
                    outputFile: outputFile, header: header,
                    pairedEnd: pe2
                )
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
                // Use best region for mate info, then emit with full supplementary handling
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

                try emitSingleEndAlignments(
                    read: read1, regions: regions1,
                    outputFile: outputFile, header: header,
                    pairedEnd: pe1
                )
                try emitSingleEndAlignments(
                    read: read2, regions: regions2,
                    outputFile: outputFile, header: header,
                    pairedEnd: pe2
                )
            }
        }
    }

    // MARK: - Supplementary/Chimeric Output

    /// Classify and emit alignments for a single read with proper primary/supplementary/secondary handling.
    /// Matches bwa-mem2's mem_reg2sam() two-pass approach.
    func emitSingleEndAlignments(
        read: ReadSequence,
        regions: [MemAlnReg],
        outputFile: borrowing HTSFile,
        header: SAMHeader,
        pairedEnd: PairedEndInfo? = nil
    ) throws {
        let scoring = options.scoring
        let rgID = options.readGroupID

        // Pass 1: Classify regions into segments
        var segments: [AlnSegment] = []
        var secondaryInfos: [(rname: String, pos: Int64, isReverse: Bool,
                              cigarString: String, nm: Int32)] = []
        var nonSecondaryCount = 0

        for (regIdx, region) in regions.enumerated() {
            guard region.score >= scoring.minOutputScore else { continue }

            let cigarInfo = generateCIGAR(read: read, region: region)
            let (rid, localPos) = index.metadata.decodePosition(cigarInfo.pos)
            let rname = rid >= 0 && rid < index.metadata.annotations.count
                ? index.metadata.annotations[Int(rid)].name : "*"

            let mapq = MappingQuality.compute(
                region: region,
                allRegions: regions,
                scoring: scoring,
                readLength: Int32(read.length)
            )

            if region.secondary >= 0 {
                // True secondary: overlaps a better region
                // Check drop ratio: skip if score < parent score * 0.5
                let parentIdx = Int(region.secondary)
                if parentIdx < regions.count && region.score < regions[parentIdx].score / 2 {
                    continue
                }
                secondaryInfos.append((
                    rname: rname,
                    pos: localPos,
                    isReverse: cigarInfo.isReverse,
                    cigarString: cigarInfo.cigarString,
                    nm: cigarInfo.nm
                ))
            } else {
                // Independent region (primary or supplementary)
                let isPrimary = nonSecondaryCount == 0
                let isSupplementary = nonSecondaryCount > 0

                segments.append(AlnSegment(
                    regionIndex: regIdx,
                    cigarInfo: cigarInfo,
                    mapq: mapq,
                    isPrimary: isPrimary,
                    isSupplementary: isSupplementary,
                    isSecondary: false,
                    rname: rname,
                    localPos: localPos,
                    rid: rid
                ))
                nonSecondaryCount += 1
            }
        }

        guard !segments.isEmpty else {
            let record = try SAMOutputBuilder.buildUnmappedRecord(
                read: read, pairedEnd: pairedEnd, readGroupID: rgID
            )
            try outputFile.write(record: record, header: header)
            return
        }

        // Cap supplementary MAPQ at primary's MAPQ, except for ALT hits
        let primaryMapq = segments[0].mapq
        for i in 1..<segments.count {
            if !regions[segments[i].regionIndex].isAlt {
                segments[i].mapq = min(segments[i].mapq, primaryMapq)
            }
        }

        // Build SA tag info from all non-secondary segments
        let saSegmentInfos: [(rname: String, pos: Int64, isReverse: Bool,
                              cigarString: String, mapq: UInt8, nm: Int32)] =
            segments.map { seg in
                (rname: seg.rname, pos: seg.localPos, isReverse: seg.cigarInfo.isReverse,
                 cigarString: seg.cigarInfo.cigarString, mapq: seg.mapq,
                 nm: seg.cigarInfo.nm)
            }

        // Build XA tag from qualifying secondaries (on primary record only)
        // Use higher XA limit when ALT secondaries are present
        let hasAltSecondary = secondaryInfos.contains { sec in
            // Check if any secondary region is ALT
            regions.contains { r in
                r.secondary >= 0 && r.isAlt
                    && index.metadata.annotations.indices.contains(Int(r.rid))
                    && index.metadata.annotations[Int(r.rid)].name == sec.rname
            }
        }
        let effectiveMaxXA = hasAltSecondary
            ? Int(scoring.maxXAHitsAlt)
            : Int(scoring.maxXAHits)
        let xaTag = SAMOutputBuilder.buildXATag(
            secondaries: secondaryInfos, maxHits: effectiveMaxXA
        )

        // Pass 2: Emit records
        for (segIdx, seg) in segments.enumerated() {
            let region = regions[seg.regionIndex]

            // SA tag: present on all non-secondary segments, excluding self
            let saTag = segments.count > 1
                ? SAMOutputBuilder.buildSATag(segments: saSegmentInfos, excludeIndex: segIdx)
                : nil

            let record = try SAMOutputBuilder.buildRecord(
                read: read,
                region: region,
                allRegions: regions,
                metadata: index.metadata,
                scoring: scoring,
                cigar: seg.cigarInfo.cigar,
                nm: seg.cigarInfo.nm,
                md: seg.cigarInfo.md,
                isPrimary: seg.isPrimary,
                isSupplementary: seg.isSupplementary,
                mapqOverride: seg.mapq,
                adjustedPos: seg.cigarInfo.pos,
                pairedEnd: pairedEnd,
                saTag: saTag,
                xaTag: seg.isPrimary ? xaTag : nil,
                readGroupID: rgID
            )
            try outputFile.write(record: record, header: header)
        }
    }

    // MARK: - Private Helpers

    /// Select rescue candidate regions: primary regions with score >= bestScore - unpairedPenalty,
    /// capped at maxMatesw. Matches bwa-mem2 lines 382-385.
    private func selectRescueCandidates(
        _ regions: [MemAlnReg], scoring: ScoringParameters
    ) -> [MemAlnReg] {
        guard let bestScore = regions.first(where: { $0.secondary < 0 })?.score else {
            return []
        }
        let threshold = bestScore - scoring.unpairedPenalty
        var candidates: [MemAlnReg] = []
        for r in regions {
            guard r.secondary < 0 && r.score >= threshold else { continue }
            candidates.append(r)
            if candidates.count >= Int(scoring.maxMatesw) { break }
        }
        return candidates
    }

    /// Align all reads in parallel and return regions indexed by read position.
    private func alignAllReads(_ reads: [ReadSequence]) async -> [[MemAlnReg]] {
        let maxConcurrency = options.scoring.numThreads
        let results = await withTaskGroup(
            of: (Int, [MemAlnReg]).self,
            returning: [(Int, [MemAlnReg])].self
        ) { group in
            var nextIdx = 0
            var collected: [(Int, [MemAlnReg])] = []

            // Seed initial batch
            while nextIdx < min(maxConcurrency, reads.count) {
                let idx = nextIdx
                let read = reads[idx]
                group.addTask { [self] in
                    (idx, self.alignRead(read))
                }
                nextIdx += 1
            }

            // As each completes, launch next
            for await result in group {
                collected.append(result)
                if nextIdx < reads.count {
                    let idx = nextIdx
                    let read = reads[idx]
                    group.addTask { [self] in
                        (idx, self.alignRead(read))
                    }
                    nextIdx += 1
                }
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
        let rgID = options.readGroupID

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
                    isPrimary: true, adjustedPos: cigar.pos, pairedEnd: mappedPE,
                    readGroupID: rgID
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: unmappedRead, pairedEnd: unmappedPE, readGroupID: rgID
                )
                try outputFile.write(record: rec2, header: header)
            } else {
                let rec1 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: unmappedRead, pairedEnd: unmappedPE, readGroupID: rgID
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildRecord(
                    read: mappedRead, region: region, allRegions: mappedRegions,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar.cigar, nm: cigar.nm, md: cigar.md,
                    isPrimary: true, adjustedPos: cigar.pos, pairedEnd: mappedPE,
                    readGroupID: rgID
                )
                try outputFile.write(record: rec2, header: header)
            }
        } catch {
            fputs("[PE] Error writing mapped/unmapped pair: \(error)\n", stderr)
        }
    }
}
