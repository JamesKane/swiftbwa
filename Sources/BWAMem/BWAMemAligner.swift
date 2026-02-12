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
    /// - Parameters:
    ///   - read: The read to align
    ///   - readId: Unique identifier for deterministic hash-based tiebreaking (matches bwa-mem2's id parameter)
    nonisolated public func alignRead(_ read: ReadSequence, readId: UInt64 = 0) -> [MemAlnReg] {
        let scoring = options.scoring

        // Phase 1: Find SMEMs
        var smems = SMEMFinder.findAllSMEMs(
            query: read.bases,
            bwt: index.bwt,
            minSeedLen: scoring.minSeedLength
        )

        guard !smems.isEmpty else { return [] }

        // Phase 1.5: Re-seeding (-y) — for long reads with high-occurrence seeds,
        // re-run SMEM finding with relaxed occurrence threshold to find additional
        // shorter, more specific seeds. Matches bwa-mem2's mem_collect_intv behavior.
        if scoring.reseedLength > 0 {
            let maxOcc = Int64(scoring.maxOccurrences)
            let hasHighOcc = smems.contains { $0.count > maxOcc }
            if hasHighOcc {
                let reseeded = SMEMFinder.findAllSMEMs(
                    query: read.bases,
                    bwt: index.bwt,
                    minSeedLen: scoring.reseedLength,
                    minIntv: maxOcc
                )
                // Merge and deduplicate: keep unique seeds by (queryBegin, queryEnd, k)
                var seen = Set<Int64>()
                for s in smems {
                    // Simple hash combining position and SA interval
                    seen.insert(Int64(s.queryBegin) << 32 | Int64(s.queryEnd) << 16 | (s.k & 0xFFFF))
                }
                for s in reseeded {
                    let key = Int64(s.queryBegin) << 32 | Int64(s.queryEnd) << 16 | (s.k & 0xFFFF)
                    if !seen.contains(key) {
                        smems.append(s)
                        seen.insert(key)
                    }
                }
                // Re-sort
                smems.sort {
                    if $0.queryBegin != $1.queryBegin { return $0.queryBegin < $1.queryBegin }
                    return $0.length > $1.length
                }
            }
        }

        // Phase 2: Chain seeds
        var chains = SeedChainer.chain(
            smems: smems,
            getSAEntry: { [index] pos in
                index.suffixArray.resolve(at: pos, bwt: index.bwt)
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
                    index.getReference(at: pos, length: length)
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

        // Phase 4.5: Sort, deduplicate, and patch overlapping regions
        RegionDedup.sortDedupPatch(
            regions: &regions,
            query: read.bases,
            getReference: { [index] pos, length in
                index.getReference(at: pos, length: length)
            },
            genomeLength: index.genomeLength,
            scoring: scoring
        )

        // Phase 5: Mark secondary alignments (ALT-aware if any ALT regions)
        let hasAlt = regions.contains { $0.isAlt }
        if hasAlt {
            ChainFilter.markSecondaryALT(
                regions: &regions, maskLevel: scoring.maskLevel, scoring: scoring,
                readId: readId
            )
        } else {
            ChainFilter.markSecondary(regions: &regions, maskLevel: scoring.maskLevel,
                                      readId: readId)
        }

        return regions
    }

    /// Dummy CIGARInfo for regions filtered before CIGAR generation.
    /// Matching bwa-mem2's approach of only running ksw_global2 for output regions.
    private static let dummyCigar = CIGARInfo(
        cigar: [], nm: 0, md: "", pos: 0,
        isReverse: false, refConsumed: 0, cigarString: "*"
    )

    /// Generate CIGARs for regions, skipping those that won't be output.
    /// Regions with score < minOutputScore or secondary regions below drop ratio
    /// are replaced with a lightweight dummy, avoiding expensive GlobalAligner calls.
    nonisolated func generateFilteredCIGARs(
        read: ReadSequence,
        regions: [MemAlnReg],
        scoringMatrix: [Int8]
    ) -> [CIGARInfo] {
        let minScore = options.scoring.minOutputScore
        return regions.enumerated().map { (idx, reg) -> CIGARInfo in
            // Always generate CIGAR for primary (idx 0) — it's always emitted
            if idx > 0 && reg.score < minScore { return Self.dummyCigar }
            if reg.secondary >= 0 {
                let pi = Int(reg.secondary)
                if pi < regions.count && reg.score < regions[pi].score / 2 {
                    return Self.dummyCigar
                }
            }
            return generateCIGAR(read: read, region: reg, scoringMatrix: scoringMatrix)
        }
    }

    /// Generate CIGAR for a single alignment region.
    nonisolated func generateCIGAR(read: ReadSequence, region: MemAlnReg, scoringMatrix: [Int8]? = nil) -> CIGARInfo {
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
            refPos: rb,
            scoringMatrix: scoringMatrix
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
        // Process reads in parallel: alignment + CIGAR generation in one task
        // (CIGAR generation via GlobalAligner is the SE bottleneck, so it must be parallel)
        let maxConcurrency = options.scoring.numThreads

        let results = await withTaskGroup(
            of: (Int, [MemAlnReg], [CIGARInfo]).self,
            returning: [(Int, [MemAlnReg], [CIGARInfo])].self
        ) { group in
            var nextIdx = 0
            var collected: [(Int, [MemAlnReg], [CIGARInfo])] = []
            collected.reserveCapacity(reads.count)

            // Seed initial batch
            while nextIdx < min(maxConcurrency, reads.count) {
                let idx = nextIdx
                let read = reads[idx]
                group.addTask { [self] in
                    let regions = self.alignRead(read, readId: UInt64(idx))
                    let mat = self.options.scoring.scoringMatrix()
                    let cigars = self.generateFilteredCIGARs(read: read, regions: regions, scoringMatrix: mat)
                    return (idx, regions, cigars)
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
                        let regions = self.alignRead(read, readId: UInt64(idx))
                        let mat = self.options.scoring.scoringMatrix()
                        let cigars = self.generateFilteredCIGARs(read: read, regions: regions, scoringMatrix: mat)
                        return (idx, regions, cigars)
                    }
                    nextIdx += 1
                }
            }

            return collected  // No sort needed: output loop uses idx-based lookup
        }

        // Place results by idx for ordered output
        var regionsByIdx = Array(repeating: [MemAlnReg](), count: reads.count)
        var cigarsByIdx = Array(repeating: [CIGARInfo](), count: reads.count)
        for (idx, regions, cigars) in results {
            regionsByIdx[idx] = regions
            cigarsByIdx[idx] = cigars
        }

        // Write output sequentially (ordered by input) — CIGARs already computed
        let rgID = options.readGroupID
        let comment = options.appendComment
        for idx in 0..<reads.count {
            let read = reads[idx]
            let regions = regionsByIdx[idx]
            let cigars = cigarsByIdx[idx]

            if regions.isEmpty {
                let record = try SAMOutputBuilder.buildUnmappedRecord(
                    read: read, readGroupID: rgID,
                    appendComment: comment
                )
                try outputFile.write(record: record, header: header)
            } else {
                try emitSingleEndAlignments(
                    read: read,
                    regions: regions,
                    outputFile: outputFile,
                    header: header,
                    cigarCache: cigars
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

        // Phase 1: Align both mates concurrently
        // PE read IDs match bwa-mem2: id = pairIndex<<1|mateFlag
        async let r1Task = alignAllReads(reads1, idBase: 0, idShift: true)
        async let r2Task = alignAllReads(reads2, idBase: 1, idShift: true)
        let allRegions1 = await r1Task
        let allRegions2 = await r2Task

        // Phase 2: Estimate insert size distribution (or use manual override)
        let dist: InsertSizeDistribution
        if let manual = options.manualInsertSize {
            dist = InsertSizeEstimator.buildManualDistribution(override: manual)
        } else {
            dist = InsertSizeEstimator.estimate(
                regions1: allRegions1,
                regions2: allRegions2,
                genomeLength: genomeLen
            )
        }

        let primaryStats = dist.stats[dist.primaryOrientation.rawValue]
        if !primaryStats.failed {
            if options.verbosity >= 3 {
                fputs("[PE] Insert size: mean=\(String(format: "%.1f", primaryStats.mean)), "
                      + "stddev=\(String(format: "%.1f", primaryStats.stddev)), "
                      + "orientation=\(dist.primaryOrientation), "
                      + "n=\(primaryStats.count)\n", stderr)
            }
        } else {
            if options.verbosity >= 2 {
                fputs("[PE] Warning: insert size estimation failed, "
                      + "treating as unpaired for scoring\n", stderr)
            }
        }

        // Phase 2.5+3: Merged rescue + CIGAR pre-computation in one parallel pass.
        // Eliminates the synchronization barrier between separate rescue and CIGAR phases.
        var mutableRegions1 = allRegions1
        var mutableRegions2 = allRegions2
        let skipRescue = (options.scoring.flag & ScoringParameters.flagNoRescue) != 0
        let doRescue = !primaryStats.failed && !skipRescue

        // Batch rescue + CIGAR tasks: process multiple pairs per task to reduce
        // scheduling overhead (500k individual tasks → ~4k batched tasks).
        let maxConcurrencyOut = options.scoring.numThreads
        let batchSize = 128
        let batchCount = (pairCount + batchSize - 1) / batchSize

        typealias PairResult = (Int, [MemAlnReg], [MemAlnReg], [CIGARInfo], [CIGARInfo], Int, Int)
        let mergedResults = await withTaskGroup(
            of: [PairResult].self,
            returning: [[PairResult]].self
        ) { group in
            let scoring = options.scoring
            // Capture arrays for task access (CoW reference counted)
            let capturedReads1 = reads1
            let capturedReads2 = reads2
            let capturedRegs1 = mutableRegions1
            let capturedRegs2 = mutableRegions2

            var nextBatch = 0
            var collected: [[PairResult]] = []
            collected.reserveCapacity(batchCount)

            func addBatchTask(
                _ group: inout TaskGroup<[PairResult]>,
                batchStart: Int
            ) {
                let batchEnd = min(batchStart + batchSize, pairCount)
                group.addTask { [self] in
                    let mat = scoring.scoringMatrix()
                    var results: [PairResult] = []
                    results.reserveCapacity(batchEnd - batchStart)
                    for idx in batchStart..<batchEnd {
                        let r1 = capturedReads1[idx]
                        let r2 = capturedReads2[idx]
                        var finalRegs1 = capturedRegs1[idx]
                        var finalRegs2 = capturedRegs2[idx]
                        var rc1 = 0, rc2 = 0
                        if doRescue {
                            let (_, rr1, rr2, c1, c2) = self.rescuePair(
                                idx: idx, regions1: finalRegs1, regions2: finalRegs2,
                                read1: r1, read2: r2, dist: dist, scoring: scoring
                            )
                            finalRegs1 = rr1; finalRegs2 = rr2; rc1 = c1; rc2 = c2
                        }
                        let cig1 = self.generateFilteredCIGARs(read: r1, regions: finalRegs1, scoringMatrix: mat)
                        let cig2 = self.generateFilteredCIGARs(read: r2, regions: finalRegs2, scoringMatrix: mat)
                        results.append((idx, finalRegs1, finalRegs2, cig1, cig2, rc1, rc2))
                    }
                    return results
                }
            }

            // Seed initial batches
            while nextBatch < min(maxConcurrencyOut, batchCount) {
                addBatchTask(&group, batchStart: nextBatch * batchSize)
                nextBatch += 1
            }

            for await batchResult in group {
                collected.append(batchResult)
                if nextBatch < batchCount {
                    addBatchTask(&group, batchStart: nextBatch * batchSize)
                    nextBatch += 1
                }
            }

            return collected  // No sort needed: unpacking uses idx-based placement
        }

        // Unpack merged results (batched)
        var allCigars = Array(repeating: ([CIGARInfo](), [CIGARInfo]()), count: pairCount)
        var rescueCount1 = 0, rescueCount2 = 0
        for batch in mergedResults {
            for (idx, regs1, regs2, cig1, cig2, rc1, rc2) in batch {
                mutableRegions1[idx] = regs1
                mutableRegions2[idx] = regs2
                allCigars[idx] = (cig1, cig2)
                rescueCount1 += rc1
                rescueCount2 += rc2
            }
        }

        if doRescue && options.verbosity >= 3 {
            fputs("[PE] Mate rescue: \(rescueCount1) read1 + \(rescueCount2) read2 rescued\n", stderr)
        }

        let skipPairing = (options.scoring.flag & ScoringParameters.flagNoPairing) != 0
        let rgID = options.readGroupID
        let comment = options.appendComment
        for i in 0..<pairCount {
            let read1 = reads1[i]
            let read2 = reads2[i]
            let regions1 = mutableRegions1[i]
            let regions2 = mutableRegions2[i]
            let (cigars1, cigars2) = allCigars[i]

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
                    read: read1, pairedEnd: pe1, readGroupID: rgID,
                    appendComment: comment
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: read2, pairedEnd: pe2, readGroupID: rgID,
                    appendComment: comment
                )
                try outputFile.write(record: rec2, header: header)
                continue
            }

            // Try to resolve best pair (skip if -P)
            let decision = skipPairing ? nil : PairedEndResolver.resolve(
                regions1: regions1,
                regions2: regions2,
                dist: dist,
                genomeLength: genomeLen,
                scoring: options.scoring
            )

            if let decision = decision {
                // Both mapped and paired
                // Promote paired regions to primary position (matches bwa-mem2's
                // mem_sam_pe swap logic that puts the paired region at index 0)
                var pairedRegions1 = regions1
                var pairedRegions2 = regions2
                var pairedCigars1 = cigars1
                var pairedCigars2 = cigars2
                promotePairedRegion(regions: &pairedRegions1, pairedIdx: decision.idx1)
                promotePairedRegion(regions: &pairedRegions2, pairedIdx: decision.idx2)
                if decision.idx1 > 0 && decision.idx1 < pairedCigars1.count {
                    pairedCigars1.swapAt(0, decision.idx1)
                }
                if decision.idx2 > 0 && decision.idx2 < pairedCigars2.count {
                    pairedCigars2.swapAt(0, decision.idx2)
                }

                let cigar1 = pairedCigars1[0]
                let cigar2 = pairedCigars2[0]

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
                    read: read1, regions: pairedRegions1,
                    outputFile: outputFile, header: header,
                    pairedEnd: pe1, cigarCache: pairedCigars1
                )

                // Emit primary read2 using emitSingleEndAlignments for full supplementary handling
                try emitSingleEndAlignments(
                    read: read2, regions: pairedRegions2,
                    outputFile: outputFile, header: header,
                    pairedEnd: pe2, cigarCache: pairedCigars2
                )
            } else if !r1Empty && r2Empty {
                // Read 1 mapped, read 2 unmapped
                writeMappedUnmappedPair(
                    mappedRead: read1, mappedRegions: regions1,
                    unmappedRead: read2,
                    mappedIsRead1: true,
                    outputFile: outputFile, header: header,
                    cigarCache: cigars1
                )
            } else if r1Empty && !r2Empty {
                // Read 1 unmapped, read 2 mapped
                writeMappedUnmappedPair(
                    mappedRead: read2, mappedRegions: regions2,
                    unmappedRead: read1,
                    mappedIsRead1: false,
                    outputFile: outputFile, header: header,
                    cigarCache: cigars2
                )
            } else {
                // Both have regions but no valid pairing (e.g., different chromosomes)
                // Use best region for mate info, then emit with full supplementary handling
                let cigar1 = cigars1[0]
                let cigar2 = cigars2[0]

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
                    pairedEnd: pe1, cigarCache: cigars1
                )
                try emitSingleEndAlignments(
                    read: read2, regions: regions2,
                    outputFile: outputFile, header: header,
                    pairedEnd: pe2, cigarCache: cigars2
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
        pairedEnd: PairedEndInfo? = nil,
        cigarCache: [CIGARInfo]? = nil
    ) throws {
        let scoring = options.scoring
        let rgID = options.readGroupID
        let refHeader = options.outputRefHeader
        let comment = options.appendComment
        let outputAll = (scoring.flag & ScoringParameters.flagAll) != 0
        let scoringMat = scoring.scoringMatrix()

        // Pass 1: Classify regions into segments
        var segments: [AlnSegment] = []
        var secondaryInfos: [(rname: String, pos: Int64, isReverse: Bool,
                              cigarString: String, nm: Int32)] = []
        var nonSecondaryCount = 0

        for (regIdx, region) in regions.enumerated() {
            guard region.score >= scoring.minOutputScore else { continue }

            let cigarInfo: CIGARInfo
            if let cache = cigarCache, regIdx < cache.count {
                cigarInfo = cache[regIdx]
            } else {
                cigarInfo = generateCIGAR(read: read, region: region, scoringMatrix: scoringMat)
            }
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
                if outputAll {
                    // -a: emit as full SAM record with 0x100 flag
                    segments.append(AlnSegment(
                        regionIndex: regIdx,
                        cigarInfo: cigarInfo,
                        mapq: mapq,
                        isPrimary: false,
                        isSupplementary: false,
                        isSecondary: true,
                        rname: rname,
                        localPos: localPos,
                        rid: rid
                    ))
                } else {
                    secondaryInfos.append((
                        rname: rname,
                        pos: localPos,
                        isReverse: cigarInfo.isReverse,
                        cigarString: cigarInfo.cigarString,
                        nm: cigarInfo.nm
                    ))
                }
            } else {
                // Independent region (primary or supplementary)
                let noMulti = (scoring.flag & ScoringParameters.flagNoMulti) != 0
                let isPrimary = nonSecondaryCount == 0
                let isSupplementary = nonSecondaryCount > 0 && !noMulti
                let isSecondary = nonSecondaryCount > 0 && noMulti

                segments.append(AlnSegment(
                    regionIndex: regIdx,
                    cigarInfo: cigarInfo,
                    mapq: mapq,
                    isPrimary: isPrimary,
                    isSupplementary: isSupplementary,
                    isSecondary: isSecondary,
                    rname: rname,
                    localPos: localPos,
                    rid: rid
                ))
                nonSecondaryCount += 1
            }
        }

        guard !segments.isEmpty else {
            let record = try SAMOutputBuilder.buildUnmappedRecord(
                read: read, pairedEnd: pairedEnd, readGroupID: rgID,
                appendComment: comment
            )
            try outputFile.write(record: record, header: header)
            return
        }

        // -5: Reorder so smallest-coordinate segment becomes primary
        if (scoring.flag & ScoringParameters.flagPrimary5) != 0 && segments.count > 1 {
            var minIdx = 0
            for i in 1..<segments.count {
                if segments[i].localPos < segments[minIdx].localPos {
                    minIdx = i
                }
            }
            if minIdx != 0 {
                segments.swapAt(0, minIdx)
                for i in 0..<segments.count {
                    segments[i] = AlnSegment(
                        regionIndex: segments[i].regionIndex,
                        cigarInfo: segments[i].cigarInfo,
                        mapq: segments[i].mapq,
                        isPrimary: i == 0,
                        isSupplementary: i > 0,
                        isSecondary: segments[i].isSecondary,
                        rname: segments[i].rname,
                        localPos: segments[i].localPos,
                        rid: segments[i].rid
                    )
                }
            }
        }

        // Cap supplementary MAPQ at primary's MAPQ, except for ALT hits
        // Skip if -q or -5 (flagKeepSuppMapq)
        if (scoring.flag & ScoringParameters.flagKeepSuppMapq) == 0 {
            let primaryMapq = segments[0].mapq
            for i in 1..<segments.count {
                if !regions[segments[i].regionIndex].isAlt {
                    segments[i].mapq = min(segments[i].mapq, primaryMapq)
                }
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
        // When -a is set, secondaries are emitted as full records, so skip XA
        let xaTag: String?
        if outputAll {
            xaTag = nil
        } else {
            let hasAltSecondary = secondaryInfos.contains { sec in
                regions.contains { r in
                    r.secondary >= 0 && r.isAlt
                        && index.metadata.annotations.indices.contains(Int(r.rid))
                        && index.metadata.annotations[Int(r.rid)].name == sec.rname
                }
            }
            let effectiveMaxXA = hasAltSecondary
                ? Int(scoring.maxXAHitsAlt)
                : Int(scoring.maxXAHits)
            xaTag = SAMOutputBuilder.buildXATag(
                secondaries: secondaryInfos, maxHits: effectiveMaxXA
            )
        }

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
                readGroupID: rgID,
                outputRefHeader: refHeader,
                appendComment: comment
            )
            try outputFile.write(record: record, header: header)
        }
    }

    // MARK: - Private Helpers

    /// Promote the paired region to index 0 for primary output.
    /// Matches bwa-mem2's mem_sam_pe() swap logic that ensures the paired region
    /// is emitted as primary. Fixes secondary references accordingly.
    private func promotePairedRegion(regions: inout [MemAlnReg], pairedIdx: Int) {
        guard pairedIdx > 0 && pairedIdx < regions.count else { return }

        // Swap the paired region to position 0
        regions.swapAt(0, pairedIdx)

        // Fix secondary references that point to the swapped indices
        for i in 0..<regions.count {
            if regions[i].secondary == Int32(pairedIdx) {
                regions[i].secondary = 0
            } else if regions[i].secondary == 0 {
                regions[i].secondary = Int32(pairedIdx)
            }
        }

        // Ensure the promoted region is marked as primary
        regions[0].secondary = -1
    }

    /// Rescue a single pair: rescue read2 from read1's templates, then read1 from
    /// (read2 + rescued). Re-marks secondaries. Called in parallel across pairs.
    nonisolated private func rescuePair(
        idx: Int,
        regions1: [MemAlnReg],
        regions2: [MemAlnReg],
        read1: ReadSequence,
        read2: ReadSequence,
        dist: InsertSizeDistribution,
        scoring: ScoringParameters
    ) -> (Int, [MemAlnReg], [MemAlnReg], Int, Int) {
        let genomeLen = index.genomeLength
        let mat = scoring.scoringMatrix()
        var regs1 = regions1
        var regs2 = regions2

        // Rescue read2 using read1's alignments as templates
        let templates1 = Self.selectRescueCandidates(regs1, scoring: scoring)
        let rescued2 = MateRescue.rescue(
            templateRegions: templates1,
            mateRead: read2,
            mateRegions: regs2,
            dist: dist,
            genomeLength: genomeLen,
            packedRef: index.packedRef,
            metadata: index.metadata,
            scoring: scoring,
            scoringMatrix: mat
        )
        regs2.append(contentsOf: rescued2)

        // Rescue read1 using read2's alignments (symmetric)
        let templates2 = Self.selectRescueCandidates(regs2, scoring: scoring)
        let rescued1 = MateRescue.rescue(
            templateRegions: templates2,
            mateRead: read1,
            mateRegions: regs1,
            dist: dist,
            genomeLength: genomeLen,
            packedRef: index.packedRef,
            metadata: index.metadata,
            scoring: scoring,
            scoringMatrix: mat
        )
        regs1.append(contentsOf: rescued1)

        // Re-mark secondaries after adding rescued regions
        let peId1 = (UInt64(idx) &<< 1) | 0
        let peId2 = (UInt64(idx) &<< 1) | 1
        if regs1.count > 1 {
            ChainFilter.markSecondary(
                regions: &regs1, maskLevel: scoring.maskLevel, readId: peId1
            )
        }
        if regs2.count > 1 {
            ChainFilter.markSecondary(
                regions: &regs2, maskLevel: scoring.maskLevel, readId: peId2
            )
        }

        return (idx, regs1, regs2, rescued1.count, rescued2.count)
    }

    /// Select rescue candidate regions: primary regions with score >= bestScore - unpairedPenalty,
    /// capped at maxMatesw. Matches bwa-mem2 lines 382-385.
    private static func selectRescueCandidates(
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
    /// - Parameters:
    ///   - reads: Array of reads to align
    ///   - idBase: Base value for read ID (0 for SE/read1, 1 for read2)
    ///   - idShift: If true, use PE scheme (idx<<1|idBase); if false, use SE scheme (idx)
    nonisolated private func alignAllReads(
        _ reads: [ReadSequence],
        idBase: UInt64 = 0,
        idShift: Bool = false
    ) async -> [[MemAlnReg]] {
        let maxConcurrency = options.scoring.numThreads
        let results = await withTaskGroup(
            of: (Int, [MemAlnReg]).self,
            returning: [(Int, [MemAlnReg])].self
        ) { group in
            var nextIdx = 0
            var collected: [(Int, [MemAlnReg])] = []
            collected.reserveCapacity(reads.count)

            // Seed initial batch
            while nextIdx < min(maxConcurrency, reads.count) {
                let idx = nextIdx
                let read = reads[idx]
                let readId = idShift ? (UInt64(idx) &<< 1) | idBase : UInt64(idx)
                group.addTask { [self] in
                    (idx, self.alignRead(read, readId: readId))
                }
                nextIdx += 1
            }

            // As each completes, launch next
            for await result in group {
                collected.append(result)
                if nextIdx < reads.count {
                    let idx = nextIdx
                    let read = reads[idx]
                    let readId = idShift ? (UInt64(idx) &<< 1) | idBase : UInt64(idx)
                    group.addTask { [self] in
                        (idx, self.alignRead(read, readId: readId))
                    }
                    nextIdx += 1
                }
            }

            return collected  // No sort needed: idx-based placement below
        }
        // Place results by idx for ordered output
        var ordered = Array(repeating: [MemAlnReg](), count: reads.count)
        for (idx, regions) in results {
            ordered[idx] = regions
        }
        return ordered
    }

    /// Write a pair where one read is mapped and the other is unmapped.
    private func writeMappedUnmappedPair(
        mappedRead: ReadSequence,
        mappedRegions: [MemAlnReg],
        unmappedRead: ReadSequence,
        mappedIsRead1: Bool,
        outputFile: borrowing HTSFile,
        header: SAMHeader,
        cigarCache: [CIGARInfo]? = nil
    ) {
        let region = mappedRegions[0]
        let cigar = cigarCache?.first ?? generateCIGAR(read: mappedRead, region: region)
        let (rid, localPos) = index.metadata.decodePosition(cigar.pos)
        let rgID = options.readGroupID
        let refHeader = options.outputRefHeader
        let comment = options.appendComment

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
                    readGroupID: rgID,
                    outputRefHeader: refHeader, appendComment: comment
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: unmappedRead, pairedEnd: unmappedPE, readGroupID: rgID,
                    appendComment: comment
                )
                try outputFile.write(record: rec2, header: header)
            } else {
                let rec1 = try SAMOutputBuilder.buildUnmappedRecord(
                    read: unmappedRead, pairedEnd: unmappedPE, readGroupID: rgID,
                    appendComment: comment
                )
                try outputFile.write(record: rec1, header: header)
                let rec2 = try SAMOutputBuilder.buildRecord(
                    read: mappedRead, region: region, allRegions: mappedRegions,
                    metadata: index.metadata, scoring: options.scoring,
                    cigar: cigar.cigar, nm: cigar.nm, md: cigar.md,
                    isPrimary: true, adjustedPos: cigar.pos, pairedEnd: mappedPE,
                    readGroupID: rgID,
                    outputRefHeader: refHeader, appendComment: comment
                )
                try outputFile.write(record: rec2, header: header)
            }
        } catch {
            if options.verbosity >= 1 {
                fputs("[PE] Error writing mapped/unmapped pair: \(error)\n", stderr)
            }
        }
    }
}
