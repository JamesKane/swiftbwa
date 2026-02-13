import BWACore
import FMIndex
import Alignment
import Htslib
#if canImport(Metal)
import MetalSW
#endif
#if canImport(Darwin)
import Darwin
import Dispatch
#elseif canImport(Glibc)
import Glibc
import Dispatch
#elseif canImport(Musl)
import Musl
import Dispatch
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

    /// Align a single read against the reference (CPU path).
    /// Calls phase1to3, CPU extension, then phase4to5.
    nonisolated public func alignRead(_ read: ReadSequence, readId: UInt64 = 0) -> [MemAlnReg] {
        let chains = alignReadPhase1to3(read)
        guard !chains.isEmpty else { return [] }

        let scoring = options.scoring
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

        return alignReadPhase4to5(read, readId: readId, regions: regions)
    }

    /// Phase 1-3: SMEM finding → chaining → filtering. Returns filtered chains.
    nonisolated func alignReadPhase1to3(_ read: ReadSequence) -> [MemChain] {
        let scoring = options.scoring

        // Phase 1: Find SMEMs
        var smems = SMEMFinder.findAllSMEMs(
            query: read.bases,
            bwt: index.bwt,
            minSeedLen: scoring.minSeedLength
        )

        guard !smems.isEmpty else { return [] }

        // Phase 1.5: Re-seeding (-y)
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
                var seen = Set<Int64>()
                for s in smems {
                    seen.insert(Int64(s.queryBegin) << 32 | Int64(s.queryEnd) << 16 | (s.k & 0xFFFF))
                }
                for s in reseeded {
                    let key = Int64(s.queryBegin) << 32 | Int64(s.queryEnd) << 16 | (s.k & 0xFFFF)
                    if !seen.contains(key) {
                        smems.append(s)
                        seen.insert(key)
                    }
                }
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

        return chains
    }

    /// Phase 1.5-3: Re-seeding (if needed) + chaining + filtering from pre-computed SMEMs.
    /// Used by GPU seeding path where phase 1 (SMEM finding) was done on GPU.
    nonisolated func alignReadPhase1_5to3(_ read: ReadSequence, gpuSMEMs: [SMEM]) -> [MemChain] {
        let scoring = options.scoring
        var smems = gpuSMEMs
        guard !smems.isEmpty else { return [] }

        // Phase 1.5: Re-seeding (-y) — runs on CPU for high-occ seeds
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
                var seen = Set<Int64>()
                for s in smems {
                    seen.insert(Int64(s.queryBegin) << 32 | Int64(s.queryEnd) << 16 | (s.k & 0xFFFF))
                }
                for s in reseeded {
                    let key = Int64(s.queryBegin) << 32 | Int64(s.queryEnd) << 16 | (s.k & 0xFFFF)
                    if !seen.contains(key) {
                        smems.append(s)
                        seen.insert(key)
                    }
                }
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

        return chains
    }

    /// Phase 4.5-5: isAlt safety net, dedup, markSecondary. Takes extension results, returns final regions.
    nonisolated func alignReadPhase4to5(
        _ read: ReadSequence, readId: UInt64, regions: [MemAlnReg]
    ) -> [MemAlnReg] {
        var regions = regions
        guard !regions.isEmpty else { return [] }
        let scoring = options.scoring

        // Safety net: ensure isAlt is set from metadata
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

        #if canImport(Metal)
        // GPU seeding: batch SMEM finding on GPU, then per-read CPU pipeline
        let gpuSeedingEngine: MetalSWEngine?
        if let engine = options.useGPU ? MetalSWEngine.shared : nil,
           engine.smemForwardPipeline != nil {
            setupBWTForGPU(engine: engine)
            gpuSeedingEngine = engine
        } else {
            gpuSeedingEngine = nil
        }
        #endif

        let results: [(Int, [MemAlnReg], [CIGARInfo])]

        #if canImport(Metal)
        if let engine = gpuSeedingEngine {
            results = await alignBatchGPUSeeded(
                reads: reads, engine: engine, maxConcurrency: maxConcurrency
            )
        } else {
            results = await alignBatchCPU(reads: reads, maxConcurrency: maxConcurrency)
        }
        #else
        results = await alignBatchCPU(reads: reads, maxConcurrency: maxConcurrency)
        #endif

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

    /// SE CPU path: per-read parallel align + CIGAR.
    private func alignBatchCPU(
        reads: [ReadSequence],
        maxConcurrency: Int
    ) async -> [(Int, [MemAlnReg], [CIGARInfo])] {
        await withTaskGroup(
            of: (Int, [MemAlnReg], [CIGARInfo]).self,
            returning: [(Int, [MemAlnReg], [CIGARInfo])].self
        ) { group in
            var nextIdx = 0
            var collected: [(Int, [MemAlnReg], [CIGARInfo])] = []
            collected.reserveCapacity(reads.count)

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

            return collected
        }
    }

    #if canImport(Metal)
    /// SE GPU-seeded path: batch SMEM finding on GPU, then per-read CPU pipeline + CIGAR.
    private func alignBatchGPUSeeded(
        reads: [ReadSequence],
        engine: MetalSWEngine,
        maxConcurrency: Int
    ) async -> [(Int, [MemAlnReg], [CIGARInfo])] {
        let scoring = options.scoring
        let gpuBatchSize = 4096
        let readCount = reads.count
        var allCollected: [(Int, [MemAlnReg], [CIGARInfo])] = []
        allCollected.reserveCapacity(readCount)

        var batchStart = 0
        while batchStart < readCount {
            let batchEnd = min(batchStart + gpuBatchSize, readCount)
            let batchReads = Array(reads[batchStart..<batchEnd])
            let batchCount = batchEnd - batchStart

            // GPU Phase 1: Forward extension for all reads in batch
            let queries = batchReads.map { $0.bases }
            let fwdResults = SMEMDispatcher.dispatchBatch(
                queries: queries,
                engine: engine
            )
            let gpuSMEMs = SMEMDispatcher.extractSMEMs(
                from: fwdResults,
                queries: queries,
                minSeedLen: scoring.minSeedLength
            )

            // Phases 1.5-5 + CIGAR: parallel per-read on CPU
            let batchResults = await withTaskGroup(
                of: (Int, [MemAlnReg], [CIGARInfo]).self,
                returning: [(Int, [MemAlnReg], [CIGARInfo])].self
            ) { group in
                var nextIdx = 0
                var collected: [(Int, [MemAlnReg], [CIGARInfo])] = []
                collected.reserveCapacity(batchCount)

                let capturedBatchStart = batchStart
                func addTask(_ group: inout TaskGroup<(Int, [MemAlnReg], [CIGARInfo])>, idx: Int) {
                    let read = batchReads[idx]
                    let smems = gpuSMEMs[idx]
                    let gi = capturedBatchStart + idx
                    group.addTask { [self] in
                        let chains = self.alignReadPhase1_5to3(read, gpuSMEMs: smems)
                        guard !chains.isEmpty else {
                            return (gi, [], [])
                        }

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

                        let finalRegions = self.alignReadPhase4to5(
                            read, readId: UInt64(gi), regions: regions
                        )
                        let mat = scoring.scoringMatrix()
                        let cigars = self.generateFilteredCIGARs(
                            read: read, regions: finalRegions, scoringMatrix: mat
                        )
                        return (gi, finalRegions, cigars)
                    }
                }

                while nextIdx < min(maxConcurrency, batchCount) {
                    addTask(&group, idx: nextIdx)
                    nextIdx += 1
                }

                for await result in group {
                    collected.append(result)
                    if nextIdx < batchCount {
                        addTask(&group, idx: nextIdx)
                        nextIdx += 1
                    }
                }

                return collected
            }

            allCollected.append(contentsOf: batchResults)
            batchStart = batchEnd
        }

        return allCollected
    }
    #endif

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
        let allRegions1: [[MemAlnReg]]
        let allRegions2: [[MemAlnReg]]
        #if canImport(Metal)
        if let engine = options.useGPU ? MetalSWEngine.shared : nil,
           engine.smemForwardPipeline != nil {
            setupBWTForGPU(engine: engine)
            async let r1Task = alignAllReadsGPUSeeded(reads1, idBase: 0, idShift: true, engine: engine)
            async let r2Task = alignAllReadsGPUSeeded(reads2, idBase: 1, idShift: true, engine: engine)
            allRegions1 = await r1Task
            allRegions2 = await r2Task
        } else {
            async let r1Task = alignAllReads(reads1, idBase: 0, idShift: true)
            async let r2Task = alignAllReads(reads2, idBase: 1, idShift: true)
            allRegions1 = await r1Task
            allRegions2 = await r2Task
        }
        #else
        do {
            async let r1Task = alignAllReads(reads1, idBase: 0, idShift: true)
            async let r2Task = alignAllReads(reads2, idBase: 1, idShift: true)
            allRegions1 = await r1Task
            allRegions2 = await r2Task
        }
        #endif

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

        // Phase 2.5+3: Rescue + CIGAR pre-computation (merged into one parallel phase).
        var mutableRegions1 = allRegions1
        var mutableRegions2 = allRegions2
        let skipRescue = (options.scoring.flag & ScoringParameters.flagNoRescue) != 0
        let doRescue = !primaryStats.failed && !skipRescue

        let maxConcurrencyOut = options.scoring.numThreads
        #if canImport(Metal)
        // Larger batches for GPU to amortize Metal command buffer overhead
        // (~700 rescue candidates per batch → efficient GPU dispatch).
        let batchSize = (options.useGPU && MetalSWEngine.shared != nil) ? 4096 : 128
        #else
        let batchSize = 128
        #endif
        let batchCount = (pairCount + batchSize - 1) / batchSize

        var rescueCount1 = 0, rescueCount2 = 0

        #if canImport(Metal)
        let gpuEngine: MetalSWEngine? = options.useGPU ? MetalSWEngine.shared : nil
        #else
        let gpuEngine: Bool? = nil
        #endif

        typealias PairResult = (Int, [MemAlnReg], [MemAlnReg], [CIGARInfo], [CIGARInfo], Int, Int)
        let mergedResults = await withTaskGroup(
            of: [PairResult].self,
            returning: [[PairResult]].self
        ) { group in
            let scoring = options.scoring
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
                    self.processRescueCIGARBatch(
                        batchStart: batchStart, batchEnd: batchEnd,
                        reads1: capturedReads1, reads2: capturedReads2,
                        regs1: capturedRegs1, regs2: capturedRegs2,
                        doRescue: doRescue, dist: dist, scoring: scoring,
                        gpuEngine: gpuEngine
                    )
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

            return collected
        }

        // Unpack merged results
        var allCigars = Array(repeating: ([CIGARInfo](), [CIGARInfo]()), count: pairCount)
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

    /// Process a batch of pairs: rescue (CPU or GPU) + CIGAR generation.
    /// When gpuEngine is available, collects all rescue candidates for the batch,
    /// dispatches 2 GPU batches (read2 then read1), applies results, then generates CIGARs.
    /// This allows GPU rescue and CPU CIGAR generation to overlap across concurrent batches.
    nonisolated private func processRescueCIGARBatch(
        batchStart: Int, batchEnd: Int,
        reads1: [ReadSequence], reads2: [ReadSequence],
        regs1: [[MemAlnReg]], regs2: [[MemAlnReg]],
        doRescue: Bool, dist: InsertSizeDistribution,
        scoring: ScoringParameters,
        gpuEngine: Any?
    ) -> [(Int, [MemAlnReg], [MemAlnReg], [CIGARInfo], [CIGARInfo], Int, Int)] {
        let mat = scoring.scoringMatrix()
        let count = batchEnd - batchStart
        var finalRegs1 = (batchStart..<batchEnd).map { regs1[$0] }
        var finalRegs2 = (batchStart..<batchEnd).map { regs2[$0] }
        var rescueCounts1 = [Int](repeating: 0, count: count)
        var rescueCounts2 = [Int](repeating: 0, count: count)

        if doRescue {
            #if canImport(Metal)
            if let engine = gpuEngine as? MetalSWEngine {
                // GPU per-batch rescue: collect candidates, dispatch, apply
                let genomeLen = index.genomeLength

                // --- Rescue read2 from read1 templates ---
                var allCands2: [MateRescue.RescueCandidate] = []
                var candRanges2 = [(start: Int, count: Int)](repeating: (0, 0), count: count)
                for li in 0..<count {
                    let gi = batchStart + li
                    let templates = Self.selectRescueCandidates(finalRegs1[li], scoring: scoring)
                    guard !templates.isEmpty else { continue }
                    let cands = MateRescue.collectCandidates(
                        templateRegions: templates,
                        mateRead: reads2[gi],
                        mateRegions: finalRegs2[li],
                        dist: dist,
                        genomeLength: genomeLen,
                        packedRef: index.packedRef,
                        metadata: index.metadata
                    )
                    if !cands.isEmpty {
                        candRanges2[li] = (start: allCands2.count, count: cands.count)
                        allCands2.append(contentsOf: cands)
                    }
                }

                if !allCands2.isEmpty {
                    let tasks = allCands2.map {
                        LocalSWTask(query: $0.mateSeq, target: $0.refBases, scoring: scoring)
                    }
                    let gpuResults = LocalSWDispatcher.dispatchBatch(tasks: tasks, engine: engine)
                    for li in 0..<count {
                        let r = candRanges2[li]
                        for j in 0..<r.count {
                            let ci = r.start + j
                            let cand = allCands2[ci]
                            let sw: (score: Int32, qb: Int32, qe: Int32, tb: Int32, te: Int32)?
                            if let g = gpuResults[ci] {
                                sw = (g.score, g.queryBegin, g.queryEnd, g.targetBegin, g.targetEnd)
                            } else if let c = LocalSWAligner.align(
                                query: cand.mateSeq, target: cand.refBases,
                                scoring: scoring, scoringMatrix: mat
                            ) {
                                sw = (c.score, c.queryBegin, c.queryEnd, c.targetBegin, c.targetEnd)
                            } else { sw = nil }
                            if let sw = sw, let reg = MateRescue.applyResult(
                                candidate: cand, score: sw.score,
                                queryBegin: sw.qb, queryEnd: sw.qe,
                                targetBegin: sw.tb, targetEnd: sw.te,
                                genomeLength: genomeLen, scoring: scoring
                            ) {
                                finalRegs2[li].append(reg)
                                rescueCounts2[li] += 1
                            }
                        }
                    }
                }

                // --- Rescue read1 from (updated) read2 templates ---
                var allCands1: [MateRescue.RescueCandidate] = []
                var candRanges1 = [(start: Int, count: Int)](repeating: (0, 0), count: count)
                for li in 0..<count {
                    let gi = batchStart + li
                    let templates = Self.selectRescueCandidates(finalRegs2[li], scoring: scoring)
                    guard !templates.isEmpty else { continue }
                    let cands = MateRescue.collectCandidates(
                        templateRegions: templates,
                        mateRead: reads1[gi],
                        mateRegions: finalRegs1[li],
                        dist: dist,
                        genomeLength: genomeLen,
                        packedRef: index.packedRef,
                        metadata: index.metadata
                    )
                    if !cands.isEmpty {
                        candRanges1[li] = (start: allCands1.count, count: cands.count)
                        allCands1.append(contentsOf: cands)
                    }
                }

                if !allCands1.isEmpty {
                    let tasks = allCands1.map {
                        LocalSWTask(query: $0.mateSeq, target: $0.refBases, scoring: scoring)
                    }
                    let gpuResults = LocalSWDispatcher.dispatchBatch(tasks: tasks, engine: engine)
                    for li in 0..<count {
                        let r = candRanges1[li]
                        for j in 0..<r.count {
                            let ci = r.start + j
                            let cand = allCands1[ci]
                            let sw: (score: Int32, qb: Int32, qe: Int32, tb: Int32, te: Int32)?
                            if let g = gpuResults[ci] {
                                sw = (g.score, g.queryBegin, g.queryEnd, g.targetBegin, g.targetEnd)
                            } else if let c = LocalSWAligner.align(
                                query: cand.mateSeq, target: cand.refBases,
                                scoring: scoring, scoringMatrix: mat
                            ) {
                                sw = (c.score, c.queryBegin, c.queryEnd, c.targetBegin, c.targetEnd)
                            } else { sw = nil }
                            if let sw = sw, let reg = MateRescue.applyResult(
                                candidate: cand, score: sw.score,
                                queryBegin: sw.qb, queryEnd: sw.qe,
                                targetBegin: sw.tb, targetEnd: sw.te,
                                genomeLength: genomeLen, scoring: scoring
                            ) {
                                finalRegs1[li].append(reg)
                                rescueCounts1[li] += 1
                            }
                        }
                    }
                }

                // Re-mark secondaries after rescue
                for li in 0..<count {
                    let gi = batchStart + li
                    let peId1 = (UInt64(gi) &<< 1) | 0
                    let peId2 = (UInt64(gi) &<< 1) | 1
                    if finalRegs1[li].count > 1 {
                        ChainFilter.markSecondary(
                            regions: &finalRegs1[li], maskLevel: scoring.maskLevel, readId: peId1
                        )
                    }
                    if finalRegs2[li].count > 1 {
                        ChainFilter.markSecondary(
                            regions: &finalRegs2[li], maskLevel: scoring.maskLevel, readId: peId2
                        )
                    }
                }
            } else {
                // CPU per-pair rescue
                for li in 0..<count {
                    let gi = batchStart + li
                    let (_, rr1, rr2, c1, c2) = rescuePair(
                        idx: gi, regions1: finalRegs1[li], regions2: finalRegs2[li],
                        read1: reads1[gi], read2: reads2[gi], dist: dist, scoring: scoring
                    )
                    finalRegs1[li] = rr1; finalRegs2[li] = rr2
                    rescueCounts1[li] = c1; rescueCounts2[li] = c2
                }
            }
            #else
            for li in 0..<count {
                let gi = batchStart + li
                let (_, rr1, rr2, c1, c2) = rescuePair(
                    idx: gi, regions1: finalRegs1[li], regions2: finalRegs2[li],
                    read1: reads1[gi], read2: reads2[gi], dist: dist, scoring: scoring
                )
                finalRegs1[li] = rr1; finalRegs2[li] = rr2
                rescueCounts1[li] = c1; rescueCounts2[li] = c2
            }
            #endif
        }

        // CIGAR generation
        var results: [(Int, [MemAlnReg], [MemAlnReg], [CIGARInfo], [CIGARInfo], Int, Int)] = []
        results.reserveCapacity(count)
        for li in 0..<count {
            let gi = batchStart + li
            let cig1 = generateFilteredCIGARs(read: reads1[gi], regions: finalRegs1[li], scoringMatrix: mat)
            let cig2 = generateFilteredCIGARs(read: reads2[gi], regions: finalRegs2[li], scoringMatrix: mat)
            results.append((gi, finalRegs1[li], finalRegs2[li], cig1, cig2,
                            rescueCounts1[li], rescueCounts2[li]))
        }
        return results
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

    #if canImport(Metal)
    /// Ensure BWT + K-mer hash buffers are set up on the GPU engine. Called once lazily.
    nonisolated private func setupBWTForGPU(engine: MetalSWEngine) {
        guard engine.bwtCheckpointBuffer == nil else { return }
        let bwt = index.bwt
        guard let baseAddr = bwt.checkpoints.baseAddress else { return }
        let byteCount = bwt.checkpoints.count * MemoryLayout<CheckpointOCC>.stride
        engine.setupBWT(
            checkpointPointer: UnsafeRawPointer(baseAddr),
            checkpointByteCount: byteCount,
            counts: bwt.count,
            sentinelIndex: bwt.sentinelIndex
        )

        // Build 12-mer hash table: pre-computed BWT intervals for all 4^12 = 16M 12-mers.
        // BFS tree expansion: level 0 has 4 entries (single bases), each level extends
        // by all 4 bases, expanding 4× until level 11 fills the full 16M entries.
        // In-place expansion processes parents in reverse order so children don't
        // overwrite unprocessed parents.
        if options.verbosity >= 3 {
            fputs("[GPU] Building 12-mer hash table...\n", stderr)
        }
        let kmerK = 12
        let tableCount = 1 << (2 * kmerK)  // 4^12 = 16,777,216
        // Each entry: (k: Int64, l: Int64, s: Int64) = 24 bytes
        let entrySize = 3 * MemoryLayout<Int64>.size  // 24
        let tableByteCount = tableCount * entrySize    // 384 MB

        let tablePtr = UnsafeMutableRawPointer.allocate(
            byteCount: tableByteCount, alignment: MemoryLayout<Int64>.alignment
        )
        let table = tablePtr.bindMemory(to: Int64.self, capacity: tableCount * 3)

        // Level 0: initialize 4 entries for single bases A, C, G, T
        for base in 0..<4 {
            let k = bwt.count(for: base)
            let l = bwt.count(for: 3 - base)
            let s = bwt.count(forNext: base) - bwt.count(for: base)
            table[base * 3]     = k
            table[base * 3 + 1] = l
            table[base * 3 + 2] = s
        }

        // Levels 1 through kmerK-1: extend each parent by all 4 bases
        for level in 1..<kmerK {
            let parentCount = 1 << (2 * level)  // 4^level entries at current level

            // Parallelize the last few levels (level 5+ have 1024+ parents)
            if parentCount >= 1024 {
                // Copy parents to temp buffer to avoid in-place race conditions:
                // children of one chunk can overlap with parent indices of another chunk
                let parentBytes = parentCount * 3 * MemoryLayout<Int64>.size
                let tempParents = UnsafeMutablePointer<Int64>.allocate(capacity: parentCount * 3)
                memcpy(tempParents, table, parentBytes)

                let numThreads = min(parentCount, options.scoring.numThreads)
                let chunkSize = (parentCount + numThreads - 1) / numThreads

                DispatchQueue.concurrentPerform(iterations: numThreads) { threadIdx in
                    let start = threadIdx * chunkSize
                    let end = min(start + chunkSize, parentCount)
                    for parentIdx in start..<end {
                        let pk = tempParents[parentIdx * 3]
                        let pl = tempParents[parentIdx * 3 + 1]
                        let ps = tempParents[parentIdx * 3 + 2]

                        for base in 0..<4 {
                            let childIdx = parentIdx * 4 + base
                            if ps <= 0 {
                                table[childIdx * 3]     = 0
                                table[childIdx * 3 + 1] = 0
                                table[childIdx * 3 + 2] = 0
                            } else {
                                let ext = BackwardSearch.backwardExt(
                                    bwt: bwt,
                                    interval: (k: pl, l: pk, s: ps),
                                    base: 3 - base
                                )
                                table[childIdx * 3]     = ext.l
                                table[childIdx * 3 + 1] = ext.k
                                table[childIdx * 3 + 2] = ext.s
                            }
                        }
                    }
                }
                tempParents.deallocate()
            } else {
                // Sequential for small levels
                for parentIdx in stride(from: parentCount - 1, through: 0, by: -1) {
                    let pk = table[parentIdx * 3]
                    let pl = table[parentIdx * 3 + 1]
                    let ps = table[parentIdx * 3 + 2]

                    for base in 0..<4 {
                        let childIdx = parentIdx * 4 + base
                        if ps <= 0 {
                            table[childIdx * 3]     = 0
                            table[childIdx * 3 + 1] = 0
                            table[childIdx * 3 + 2] = 0
                        } else {
                            let ext = BackwardSearch.backwardExt(
                                bwt: bwt,
                                interval: (k: pl, l: pk, s: ps),
                                base: 3 - base
                            )
                            table[childIdx * 3]     = ext.l
                            table[childIdx * 3 + 1] = ext.k
                            table[childIdx * 3 + 2] = ext.s
                        }
                    }
                }
            }
        }

        engine.setupKmerHash(pointer: UnsafeRawPointer(tablePtr), byteCount: tableByteCount)
        tablePtr.deallocate()

        if options.verbosity >= 3 {
            fputs("[GPU] 12-mer hash table ready (\(tableByteCount / 1_048_576) MB)\n", stderr)
        }
    }

    /// GPU-seeded alignment: batch SMEM finding on GPU, then per-read CPU pipeline.
    nonisolated private func alignAllReadsGPUSeeded(
        _ reads: [ReadSequence],
        idBase: UInt64, idShift: Bool,
        engine: MetalSWEngine
    ) async -> [[MemAlnReg]] {
        let scoring = options.scoring
        let maxConcurrency = scoring.numThreads
        let gpuBatchSize = 4096
        let readCount = reads.count
        var allResults = Array(repeating: [MemAlnReg](), count: readCount)

        var batchStart = 0
        while batchStart < readCount {
            let batchEnd = min(batchStart + gpuBatchSize, readCount)
            let batchReads = Array(reads[batchStart..<batchEnd])
            let batchCount = batchEnd - batchStart

            // GPU Phase 1: Forward extension for all reads in batch
            let queries = batchReads.map { $0.bases }
            let fwdResults = SMEMDispatcher.dispatchBatch(
                queries: queries,
                engine: engine
            )
            let gpuSMEMs = SMEMDispatcher.extractSMEMs(
                from: fwdResults,
                queries: queries,
                minSeedLen: scoring.minSeedLength
            )

            // Phases 1.5-5: parallel per-read on CPU
            let capturedBatchStart = batchStart
            let batchResults = await withTaskGroup(
                of: (Int, [MemAlnReg]).self,
                returning: [(Int, [MemAlnReg])].self
            ) { group in
                var nextIdx = 0
                var collected: [(Int, [MemAlnReg])] = []
                collected.reserveCapacity(batchCount)

                while nextIdx < min(maxConcurrency, batchCount) {
                    let idx = nextIdx
                    let read = batchReads[idx]
                    let smems = gpuSMEMs[idx]
                    let gi = capturedBatchStart + idx
                    let readId = idShift ? (UInt64(gi) &<< 1) | idBase : UInt64(gi)
                    group.addTask { [self] in
                        let chains = self.alignReadPhase1_5to3(read, gpuSMEMs: smems)
                        guard !chains.isEmpty else { return (idx, []) }

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

                        return (idx, self.alignReadPhase4to5(read, readId: readId, regions: regions))
                    }
                    nextIdx += 1
                }

                for await result in group {
                    collected.append(result)
                    if nextIdx < batchCount {
                        let idx = nextIdx
                        let read = batchReads[idx]
                        let smems = gpuSMEMs[idx]
                        let gi = capturedBatchStart + idx
                        let readId = idShift ? (UInt64(gi) &<< 1) | idBase : UInt64(gi)
                        group.addTask { [self] in
                            let chains = self.alignReadPhase1_5to3(read, gpuSMEMs: smems)
                            guard !chains.isEmpty else { return (idx, []) }

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

                            return (idx, self.alignReadPhase4to5(read, readId: readId, regions: regions))
                        }
                        nextIdx += 1
                    }
                }

                return collected
            }

            for (idx, regions) in batchResults {
                allResults[capturedBatchStart + idx] = regions
            }

            batchStart = batchEnd
        }

        return allResults
    }
    #endif

    /// Align all reads in parallel and return regions indexed by read position.
    /// Uses CPU per-read parallel alignment. GPU extension is not used here because
    /// the batching overhead (sequential plan collection, synchronization barriers)
    /// exceeds the per-seed GPU computation savings — extension takes <2s fully
    /// parallelized on CPU vs ~12s with GPU batching.
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

    #if canImport(Metal)
    /// GPU-accelerated extension: process reads in batches of 4096.
    /// Phase 1-3 run parallel per-read, then extension is batched to GPU,
    /// then phase 4.5-5 run per-read.
    nonisolated private func alignAllReadsGPU(
        _ reads: [ReadSequence],
        idBase: UInt64, idShift: Bool,
        engine: MetalSWEngine
    ) async -> [[MemAlnReg]] {
        let scoring = options.scoring
        let maxConcurrency = scoring.numThreads
        let gpuBatchSize = 4096
        let readCount = reads.count
        var allResults = Array(repeating: [MemAlnReg](), count: readCount)

        // Process in batches of gpuBatchSize reads
        var batchStart = 0
        while batchStart < readCount {
            let batchEnd = min(batchStart + gpuBatchSize, readCount)
            let batchReads = Array(reads[batchStart..<batchEnd])
            let batchCount = batchEnd - batchStart

            // Phase 1-3: parallel SMEM → chain → filter
            let chainResults = await withTaskGroup(
                of: (Int, [MemChain]).self,
                returning: [(Int, [MemChain])].self
            ) { group in
                var nextIdx = 0
                var collected: [(Int, [MemChain])] = []
                collected.reserveCapacity(batchCount)

                while nextIdx < min(maxConcurrency, batchCount) {
                    let idx = nextIdx
                    let read = batchReads[idx]
                    group.addTask { [self] in
                        (idx, self.alignReadPhase1to3(read))
                    }
                    nextIdx += 1
                }
                for await result in group {
                    collected.append(result)
                    if nextIdx < batchCount {
                        let idx = nextIdx
                        let read = batchReads[idx]
                        group.addTask { [self] in
                            (idx, self.alignReadPhase1to3(read))
                        }
                        nextIdx += 1
                    }
                }
                return collected
            }

            // Place by idx
            var batchChains = Array(repeating: [MemChain](), count: batchCount)
            for (idx, chains) in chainResults {
                batchChains[idx] = chains
            }

            // Collect extension plans across all reads in this batch
            var allPlans: [SeedExtensionPlan] = []
            var planRanges = [(start: Int, count: Int)](repeating: (0, 0), count: batchCount)

            for i in 0..<batchCount {
                let chains = batchChains[i]
                guard !chains.isEmpty else { continue }
                let plans = ExtensionAligner.collectExtensionPlans(
                    chains: chains,
                    read: batchReads[i],
                    getReference: { [index] pos, length in
                        index.getReference(at: pos, length: length)
                    },
                    scoring: scoring
                )
                if !plans.isEmpty {
                    planRanges[i] = (start: allPlans.count, count: plans.count)
                    allPlans.append(contentsOf: plans)
                }
            }

            let totalPlans = allPlans.count
            if totalPlans == 0 {
                for i in 0..<batchCount {
                    let gi = batchStart + i
                    let readId = idShift ? (UInt64(gi) &<< 1) | idBase : UInt64(gi)
                    allResults[gi] = alignReadPhase4to5(
                        batchReads[i], readId: readId, regions: []
                    )
                }
                batchStart = batchEnd
                continue
            }

            // --- GPU dispatch: LEFT extensions ---
            var leftTasks: [BandedSWTask] = []
            leftTasks.reserveCapacity(totalPlans)
            for plan in allPlans {
                if !plan.leftQuery.isEmpty && !plan.leftTarget.isEmpty {
                    leftTasks.append(BandedSWTask(
                        query: plan.leftQuery, target: plan.leftTarget,
                        h0: plan.seed.score, w: scoring.bandWidth, scoring: scoring
                    ))
                } else {
                    leftTasks.append(BandedSWTask(
                        query: [], target: [], h0: 0, w: scoring.bandWidth, scoring: scoring
                    ))
                }
            }

            var leftResults: [SWResult] = dispatchBandedSWBatch(
                tasks: leftTasks, engine: engine, scoring: scoring
            )
            for i in 0..<totalPlans {
                if allPlans[i].leftQuery.isEmpty || allPlans[i].leftTarget.isEmpty {
                    leftResults[i] = SWResult()
                }
            }

            // Compute right h0 for each plan from left results
            var rightH0s: [Int32] = []
            rightH0s.reserveCapacity(totalPlans)
            for i in 0..<totalPlans {
                let plan = allPlans[i]
                if !plan.leftQuery.isEmpty && !plan.leftTarget.isEmpty {
                    rightH0s.append(leftResults[i].score)
                } else {
                    rightH0s.append(plan.seed.score)
                }
            }

            // --- GPU dispatch: RIGHT extensions ---
            var rightTasks: [BandedSWTask] = []
            rightTasks.reserveCapacity(totalPlans)
            for i in 0..<totalPlans {
                let plan = allPlans[i]
                if !plan.rightQuery.isEmpty && !plan.rightTarget.isEmpty {
                    rightTasks.append(BandedSWTask(
                        query: plan.rightQuery, target: plan.rightTarget,
                        h0: rightH0s[i], w: scoring.bandWidth, scoring: scoring
                    ))
                } else {
                    rightTasks.append(BandedSWTask(
                        query: [], target: [], h0: 0, w: scoring.bandWidth, scoring: scoring
                    ))
                }
            }

            var rightResults: [SWResult] = dispatchBandedSWBatch(
                tasks: rightTasks, engine: engine, scoring: scoring
            )
            for i in 0..<totalPlans {
                if allPlans[i].rightQuery.isEmpty || allPlans[i].rightTarget.isEmpty {
                    rightResults[i] = SWResult()
                }
            }

            // Per-read: assembleRegions → phase4to5
            for i in 0..<batchCount {
                let gi = batchStart + i
                let readId = idShift ? (UInt64(gi) &<< 1) | idBase : UInt64(gi)
                let range = planRanges[i]

                if range.count == 0 {
                    allResults[gi] = alignReadPhase4to5(
                        batchReads[i], readId: readId, regions: []
                    )
                    continue
                }

                let readPlans = Array(allPlans[range.start..<(range.start + range.count)])
                let readLeftResults = Array(leftResults[range.start..<(range.start + range.count)])
                let readRightResults = Array(rightResults[range.start..<(range.start + range.count)])

                let regions = ExtensionAligner.assembleRegions(
                    plans: readPlans,
                    leftResults: readLeftResults,
                    rightResults: readRightResults,
                    scoring: scoring
                )

                allResults[gi] = alignReadPhase4to5(
                    batchReads[i], readId: readId, regions: regions
                )
            }

            batchStart = batchEnd
        }

        return allResults
    }

    /// Dispatch a batch of banded SW tasks to GPU with 8→16 overflow handling.
    /// Returns one SWResult per task. Handles empty-task passthrough.
    nonisolated private func dispatchBandedSWBatch(
        tasks: [BandedSWTask],
        engine: MetalSWEngine,
        scoring: ScoringParameters
    ) -> [SWResult] {
        guard !tasks.isEmpty else { return [] }

        // Filter out empty tasks for GPU dispatch (they'd produce garbage)
        var realIndices: [Int] = []
        var realTasks: [BandedSWTask] = []
        for (i, task) in tasks.enumerated() {
            if !task.query.isEmpty && !task.target.isEmpty {
                realIndices.append(i)
                realTasks.append(task)
            }
        }

        // Initialize all results to zero
        var results = [SWResult](repeating: SWResult(), count: tasks.count)

        guard !realTasks.isEmpty else { return results }

        // Skip GPU for very small batches
        if realTasks.count < MetalSWEngine.minBatchSize {
            let mat = scoring.scoringMatrix()
            for (idx, task) in zip(realIndices, realTasks) {
                results[idx] = task.query.withUnsafeBufferPointer { qBuf in
                    task.target.withUnsafeBufferPointer { tBuf in
                        if let r = BandedSW8.align(
                            query: qBuf, target: tBuf, scoring: scoring,
                            w: task.w, h0: task.h0, scoringMatrix: mat
                        ) {
                            return r
                        }
                        return BandedSW16.align(
                            query: qBuf, target: tBuf, scoring: scoring,
                            w: task.w, h0: task.h0, scoringMatrix: mat
                        )
                    }
                }
            }
            return results
        }

        // 8-bit GPU dispatch
        let gpu8Results = BandedSWDispatcher.dispatchBatch8(tasks: realTasks, engine: engine)

        // Collect overflows for 16-bit fallback
        var overflowIndices: [Int] = []  // indices into realTasks
        for (i, r) in gpu8Results.enumerated() {
            if let r = r {
                results[realIndices[i]] = r
            } else {
                overflowIndices.append(i)
            }
        }

        // Handle overflows
        if !overflowIndices.isEmpty {
            let overflowTasks = overflowIndices.map { realTasks[$0] }
            let gpu16Results = BandedSWDispatcher.dispatchBatch16(
                tasks: overflowTasks, engine: engine
            )
            if gpu16Results.count == overflowTasks.count {
                for (j, oi) in overflowIndices.enumerated() {
                    results[realIndices[oi]] = gpu16Results[j]
                }
            } else {
                // CPU fallback for 16-bit failures
                let mat = scoring.scoringMatrix()
                for oi in overflowIndices {
                    let task = realTasks[oi]
                    results[realIndices[oi]] = task.query.withUnsafeBufferPointer { qBuf in
                        task.target.withUnsafeBufferPointer { tBuf in
                            BandedSW16.align(
                                query: qBuf, target: tBuf, scoring: scoring,
                                w: task.w, h0: task.h0, scoringMatrix: mat
                            )
                        }
                    }
                }
            }
        }

        return results
    }
    #endif

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
