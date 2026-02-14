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

/// Pre-computed output specification for one SAM record.
/// Built during parallel output preparation, consumed during sequential emit.
struct RecordSpec: Sendable {
    let readIsRead1: Bool
    let regionIndex: Int       // -1 for unmapped
    let cigarInfo: CIGARInfo?  // nil for unmapped
    var mapq: UInt8
    let isPrimary: Bool
    let isSupplementary: Bool
    let isSecondary: Bool
    let pairedEnd: PairedEndInfo?
    let saTag: String?
    let xaTag: String?
}

/// Complete pre-computed output plan for one read pair.
struct PairOutputPlan: Sendable {
    let regions1: [MemAlnReg]
    let regions2: [MemAlnReg]
    let records: [RecordSpec]
}

/// Per-sub-batch alignment result yielded by streaming Phase 1 functions.
struct CompletedBatch: Sendable {
    let mate: Int          // 1 or 2
    let batchIndex: Int    // sub-batch ordinal (0, 1, 2, ...)
    let batchStart: Int    // global read index where this batch starts
    let reads: [ReadSequence]
    let regions: [[MemAlnReg]]
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

        // Phase 1.25: Internal seeding — split long low-occ SMEMs to find alternative loci
        internalReseed(smems: &smems, read: read, scoring: scoring)

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

    /// Internal seeding: for long SMEMs with low occurrence, search from the
    /// midpoint with higher min_intv to find shorter seeds at alternative loci.
    /// Matches bwa-mem2's split_factor/split_width logic (bwamem.cpp:695-753).
    nonisolated func internalReseed(
        smems: inout [SMEM], read: ReadSequence, scoring: ScoringParameters
    ) {
        let splitLen = Int(Float(scoring.minSeedLength) * scoring.seedSplitRatio + 0.499)
        var newSMEMs: [SMEM] = []

        for smem in smems {
            let len = Int(smem.queryEnd - smem.queryBegin)
            if len < splitLen || smem.count > Int64(scoring.splitWidth) { continue }

            let midpoint = Int(smem.queryBegin + smem.queryEnd) / 2
            let (internalSMEMs, _) = SMEMFinder.findSMEMsAtPosition(
                query: read.bases, bwt: index.bwt,
                startPos: midpoint, minSeedLen: scoring.minSeedLength,
                minIntv: smem.count + 1
            )
            newSMEMs.append(contentsOf: internalSMEMs)
        }

        guard !newSMEMs.isEmpty else { return }
        smems.append(contentsOf: newSMEMs)
        smems.sort {
            if $0.queryBegin != $1.queryBegin { return $0.queryBegin < $1.queryBegin }
            return $0.length > $1.length
        }
    }

    /// Phase 1.5-3: Re-seeding (if needed) + chaining + filtering from pre-computed SMEMs.
    /// Used by GPU seeding path where phase 1 (SMEM finding) was done on GPU.
    nonisolated func alignReadPhase1_5to3(_ read: ReadSequence, gpuSMEMs: [SMEM]) -> [MemChain] {
        let scoring = options.scoring
        var smems = gpuSMEMs
        guard !smems.isEmpty else { return [] }

        // Internal seeding for GPU-computed SMEMs
        internalReseed(smems: &smems, read: read, scoring: scoring)

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
        let scoring = options.scoring
        let minScore = scoring.minOutputScore
        let outputAll = (scoring.flag & ScoringParameters.flagAll) != 0

        // When !outputAll, secondaries only appear in the XA tag, which is omitted
        // entirely when secondary count > maxXAHits. Pre-count to avoid generating
        // CIGARs that would be thrown away.
        let skipAllSecondaries: Bool
        if outputAll {
            skipAllSecondaries = false
        } else {
            var secondaryCount = 0
            var hasAltSecondary = false
            for (idx, reg) in regions.enumerated() {
                if idx > 0 && reg.score < minScore { continue }
                if reg.secondary >= 0 {
                    let pi = Int(reg.secondary)
                    if pi < regions.count && reg.score < regions[pi].score / 2 { continue }
                    secondaryCount += 1
                    if reg.isAlt { hasAltSecondary = true }
                }
            }
            let maxXA = hasAltSecondary ? Int(scoring.maxXAHitsAlt) : Int(scoring.maxXAHits)
            skipAllSecondaries = secondaryCount > maxXA
        }

        return regions.enumerated().map { (idx, reg) -> CIGARInfo in
            // Always generate CIGAR for primary (idx 0) — it's always emitted
            if idx > 0 && reg.score < minScore { return Self.dummyCigar }
            if reg.secondary >= 0 {
                let pi = Int(reg.secondary)
                if pi < regions.count && reg.score < regions[pi].score / 2 {
                    return Self.dummyCigar
                }
                if idx > 0 && skipAllSecondaries { return Self.dummyCigar }
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

    /// Classify alignments for one read into RecordSpecs for later emission.
    /// Extracted from emitSingleEndAlignments for parallel output preparation.
    nonisolated func classifyAlignments(
        read: ReadSequence,
        regions: [MemAlnReg],
        readIsRead1: Bool,
        pairedEnd: PairedEndInfo?,
        cigarCache: [CIGARInfo],
        scoringMat: [Int8]
    ) -> [RecordSpec] {
        let scoring = options.scoring
        let outputAll = (scoring.flag & ScoringParameters.flagAll) != 0

        guard !regions.isEmpty else {
            return [RecordSpec(
                readIsRead1: readIsRead1, regionIndex: -1, cigarInfo: nil,
                mapq: 0, isPrimary: false, isSupplementary: false,
                isSecondary: false, pairedEnd: pairedEnd, saTag: nil, xaTag: nil
            )]
        }

        // Pass 1: Classify regions into segments
        var segments: [AlnSegment] = []
        var secondaryInfos: [(rname: String, pos: Int64, isReverse: Bool,
                              cigarString: String, nm: Int32)] = []
        var nonSecondaryCount = 0

        for (regIdx, region) in regions.enumerated() {
            guard region.score >= scoring.minOutputScore else { continue }

            let cigarInfo: CIGARInfo
            if regIdx < cigarCache.count {
                cigarInfo = cigarCache[regIdx]
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
                let parentIdx = Int(region.secondary)
                if parentIdx < regions.count && region.score < regions[parentIdx].score / 2 {
                    continue
                }
                if outputAll {
                    segments.append(AlnSegment(
                        regionIndex: regIdx, cigarInfo: cigarInfo, mapq: mapq,
                        isPrimary: false, isSupplementary: false, isSecondary: true,
                        rname: rname, localPos: localPos, rid: rid
                    ))
                } else {
                    secondaryInfos.append((
                        rname: rname, pos: localPos, isReverse: cigarInfo.isReverse,
                        cigarString: cigarInfo.cigarString, nm: cigarInfo.nm
                    ))
                }
            } else {
                let noMulti = (scoring.flag & ScoringParameters.flagNoMulti) != 0
                let isPrimary = nonSecondaryCount == 0
                let isSupplementary = nonSecondaryCount > 0 && !noMulti
                let isSecondary = nonSecondaryCount > 0 && noMulti

                segments.append(AlnSegment(
                    regionIndex: regIdx, cigarInfo: cigarInfo, mapq: mapq,
                    isPrimary: isPrimary, isSupplementary: isSupplementary,
                    isSecondary: isSecondary, rname: rname, localPos: localPos, rid: rid
                ))
                nonSecondaryCount += 1
            }
        }

        guard !segments.isEmpty else {
            return [RecordSpec(
                readIsRead1: readIsRead1, regionIndex: -1, cigarInfo: nil,
                mapq: 0, isPrimary: false, isSupplementary: false,
                isSecondary: false, pairedEnd: pairedEnd, saTag: nil, xaTag: nil
            )]
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

        // Cap supplementary MAPQ at primary's MAPQ
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

        // Build XA tag from qualifying secondaries
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

        // Build RecordSpecs
        var specs: [RecordSpec] = []
        specs.reserveCapacity(segments.count)
        for (segIdx, seg) in segments.enumerated() {
            let saTag = segments.count > 1
                ? SAMOutputBuilder.buildSATag(segments: saSegmentInfos, excludeIndex: segIdx)
                : nil
            specs.append(RecordSpec(
                readIsRead1: readIsRead1,
                regionIndex: seg.regionIndex,
                cigarInfo: seg.cigarInfo,
                mapq: seg.mapq,
                isPrimary: seg.isPrimary,
                isSupplementary: seg.isSupplementary,
                isSecondary: seg.isSecondary,
                pairedEnd: pairedEnd,
                saTag: saTag,
                xaTag: seg.isPrimary ? xaTag : nil
            ))
        }
        return specs
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

        let results: [(Int, [MemAlnReg], [RecordSpec])]

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
        var specsByIdx = [[RecordSpec]?](repeating: nil, count: reads.count)
        for (idx, regions, specs) in results {
            regionsByIdx[idx] = regions
            specsByIdx[idx] = specs
        }

        let rgID = options.readGroupID
        let comment = options.appendComment
        let refHeader = options.outputRefHeader

        // Sequential record building + writing
        for idx in 0..<reads.count {
            let read = reads[idx]
            let regions = regionsByIdx[idx]
            for spec in specsByIdx[idx]! {
                if spec.regionIndex < 0 {
                    let record = try SAMOutputBuilder.buildUnmappedRecord(
                        read: read, readGroupID: rgID, appendComment: comment
                    )
                    try outputFile.write(record: record, header: header)
                } else {
                    let cigar = spec.cigarInfo!
                    let record = try SAMOutputBuilder.buildRecord(
                        read: read, region: regions[spec.regionIndex],
                        allRegions: regions, metadata: index.metadata,
                        scoring: options.scoring, cigar: cigar.cigar,
                        nm: cigar.nm, md: cigar.md,
                        isPrimary: spec.isPrimary,
                        isSupplementary: spec.isSupplementary,
                        mapqOverride: spec.mapq, adjustedPos: cigar.pos,
                        pairedEnd: nil, saTag: spec.saTag,
                        xaTag: spec.xaTag, readGroupID: rgID,
                        outputRefHeader: refHeader, appendComment: comment
                    )
                    try outputFile.write(record: record, header: header)
                }
            }
        }
    }

    /// SE CPU path: per-read parallel align + CIGAR + classify.
    private func alignBatchCPU(
        reads: [ReadSequence],
        maxConcurrency: Int
    ) async -> [(Int, [MemAlnReg], [RecordSpec])] {
        await withTaskGroup(
            of: (Int, [MemAlnReg], [RecordSpec]).self,
            returning: [(Int, [MemAlnReg], [RecordSpec])].self
        ) { group in
            var nextIdx = 0
            var collected: [(Int, [MemAlnReg], [RecordSpec])] = []
            collected.reserveCapacity(reads.count)

            while nextIdx < min(maxConcurrency, reads.count) {
                let idx = nextIdx
                let read = reads[idx]
                group.addTask { [self] in
                    let regions = self.alignRead(read, readId: UInt64(idx))
                    let mat = self.options.scoring.scoringMatrix()
                    let cigars = self.generateFilteredCIGARs(read: read, regions: regions, scoringMatrix: mat)
                    let specs = self.classifyAlignments(
                        read: read, regions: regions, readIsRead1: true,
                        pairedEnd: nil, cigarCache: cigars, scoringMat: mat
                    )
                    return (idx, regions, specs)
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
                        let specs = self.classifyAlignments(
                            read: read, regions: regions, readIsRead1: true,
                            pairedEnd: nil, cigarCache: cigars, scoringMat: mat
                        )
                        return (idx, regions, specs)
                    }
                    nextIdx += 1
                }
            }

            return collected
        }
    }

    #if canImport(Metal)
    /// SE GPU-seeded path: batch SMEM finding on GPU, then per-read CPU pipeline + CIGAR + classify.
    private func alignBatchGPUSeeded(
        reads: [ReadSequence],
        engine: MetalSWEngine,
        maxConcurrency: Int
    ) async -> [(Int, [MemAlnReg], [RecordSpec])] {
        let scoring = options.scoring
        let gpuBatchSize = 4096
        let readCount = reads.count
        var allCollected: [(Int, [MemAlnReg], [RecordSpec])] = []
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

            // Phases 1.5-5 + CIGAR + classify: parallel per-read on CPU
            let batchResults = await withTaskGroup(
                of: (Int, [MemAlnReg], [RecordSpec]).self,
                returning: [(Int, [MemAlnReg], [RecordSpec])].self
            ) { group in
                var nextIdx = 0
                var collected: [(Int, [MemAlnReg], [RecordSpec])] = []
                collected.reserveCapacity(batchCount)

                let capturedBatchStart = batchStart
                func addTask(_ group: inout TaskGroup<(Int, [MemAlnReg], [RecordSpec])>, idx: Int) {
                    let read = batchReads[idx]
                    let smems = gpuSMEMs[idx]
                    let gi = capturedBatchStart + idx
                    group.addTask { [self] in
                        let chains = self.alignReadPhase1_5to3(read, gpuSMEMs: smems)
                        guard !chains.isEmpty else {
                            let specs = self.classifyAlignments(
                                read: read, regions: [], readIsRead1: true,
                                pairedEnd: nil, cigarCache: [], scoringMat: []
                            )
                            return (gi, [], specs)
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
                        let specs = self.classifyAlignments(
                            read: read, regions: finalRegions, readIsRead1: true,
                            pairedEnd: nil, cigarCache: cigars, scoringMat: mat
                        )
                        return (gi, finalRegions, specs)
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

        // Phase 1 streams alignment results per-sub-batch via AsyncStream,
        // allowing Phase 3 to begin as soon as both mates for a sub-batch are ready.
        let (stream, continuation) = AsyncStream.makeStream(
            of: CompletedBatch.self, bufferingPolicy: .unbounded
        )

        #if canImport(Metal)
        let gpuEngine: MetalSWEngine?
        if let engine = options.useGPU ? MetalSWEngine.shared : nil,
           engine.smemForwardPipeline != nil {
            setupBWTForGPU(engine: engine)
            gpuEngine = engine
        } else {
            gpuEngine = nil
        }
        #else
        let gpuEngine: Bool? = nil
        #endif

        // Launch Phase 1 producers (both mates concurrently, yields sub-batches)
        async let phase1Done: Void = runPhase1Streaming(
            reads1: reads1, reads2: reads2,
            gpuEngine: gpuEngine,
            continuation: continuation
        )

        // Phase 2+3: Consumer pipeline — pair sub-batches, estimate insert size, process
        let skipRescue = (options.scoring.flag & ScoringParameters.flagNoRescue) != 0
        let skipPairing = (options.scoring.flag & ScoringParameters.flagNoPairing) != 0
        let rgID = options.readGroupID
        let comment = options.appendComment
        let refHeader = options.outputRefHeader
        let scoringMat = options.scoring.scoringMatrix()
        let maxConcurrencyOut = options.scoring.numThreads

        #if canImport(Metal)
        let p3BatchSize = gpuEngine != nil ? 4096 : 128
        #else
        let p3BatchSize = 128
        #endif

        var rescueCount1 = 0, rescueCount2 = 0
        var plans = [PairOutputPlan?](repeating: nil, count: pairCount)
        var nextWrite = 0

        // Pairing state: buffer mate1/mate2 sub-batch results until both arrive
        var pending1: [Int: CompletedBatch] = [:]
        var pending2: [Int: CompletedBatch] = [:]
        var insertDist: InsertSizeDistribution? = nil

        typealias FullBatchResult = [(Int, PairOutputPlan, Int, Int)]
        try await withThrowingTaskGroup(of: FullBatchResult.self) { phase3Group in
            let scoring = options.scoring
            var activeTasks = 0

            // Inline helper: process a Phase 3 result and write consecutive output
            func handlePhase3Result(_ result: FullBatchResult) throws {
                for (idx, plan, rc1, rc2) in result {
                    plans[idx] = plan
                    rescueCount1 += rc1
                    rescueCount2 += rc2
                }
                while nextWrite < pairCount, let plan = plans[nextWrite] {
                    for spec in plan.records {
                        let read = spec.readIsRead1 ? reads1[nextWrite] : reads2[nextWrite]
                        if spec.regionIndex < 0 {
                            let record = try SAMOutputBuilder.buildUnmappedRecord(
                                read: read, pairedEnd: spec.pairedEnd,
                                readGroupID: rgID, appendComment: comment
                            )
                            try outputFile.write(record: record, header: header)
                        } else {
                            let regions = spec.readIsRead1 ? plan.regions1 : plan.regions2
                            let cigar = spec.cigarInfo!
                            let record = try SAMOutputBuilder.buildRecord(
                                read: read, region: regions[spec.regionIndex],
                                allRegions: regions, metadata: index.metadata,
                                scoring: options.scoring, cigar: cigar.cigar,
                                nm: cigar.nm, md: cigar.md,
                                isPrimary: spec.isPrimary,
                                isSupplementary: spec.isSupplementary,
                                mapqOverride: spec.mapq, adjustedPos: cigar.pos,
                                pairedEnd: spec.pairedEnd, saTag: spec.saTag,
                                xaTag: spec.xaTag, readGroupID: rgID,
                                outputRefHeader: refHeader, appendComment: comment
                            )
                            try outputFile.write(record: record, header: header)
                        }
                    }
                    plans[nextWrite] = nil
                    nextWrite += 1
                }
            }

            for await batch in stream {
                // Buffer arriving sub-batch by mate
                if batch.mate == 1 { pending1[batch.batchIndex] = batch }
                else { pending2[batch.batchIndex] = batch }

                // Check if both mates are ready for this sub-batch index
                guard let m1 = pending1[batch.batchIndex],
                      let m2 = pending2[batch.batchIndex] else { continue }
                pending1.removeValue(forKey: batch.batchIndex)
                pending2.removeValue(forKey: batch.batchIndex)

                // Phase 2: Estimate insert size from first ready sub-batch pair
                if insertDist == nil {
                    if let manual = options.manualInsertSize {
                        insertDist = InsertSizeEstimator.buildManualDistribution(override: manual)
                    } else {
                        insertDist = InsertSizeEstimator.estimate(
                            regions1: m1.regions, regions2: m2.regions,
                            genomeLength: genomeLen
                        )
                    }
                    let dist = insertDist!
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
                }

                let dist = insertDist!
                let doRescue = !dist.stats[dist.primaryOrientation.rawValue].failed && !skipRescue

                // Launch Phase 3 sub-batches for this paired sub-batch
                let batchPairCount = m1.regions.count
                for offset in stride(from: 0, to: batchPairCount, by: p3BatchSize) {
                    let end = min(offset + p3BatchSize, batchPairCount)
                    let globalOff = m1.batchStart + offset
                    let subReads1 = Array(m1.reads[offset..<end])
                    let subReads2 = Array(m2.reads[offset..<end])
                    let subRegs1 = Array(m1.regions[offset..<end])
                    let subRegs2 = Array(m2.regions[offset..<end])

                    // Bounded concurrency: drain a completed Phase 3 task before adding more
                    while activeTasks >= maxConcurrencyOut {
                        if let result = try await phase3Group.next() {
                            activeTasks -= 1
                            try handlePhase3Result(result)
                        }
                    }

                    phase3Group.addTask { [self] in
                        await self.processFullBatch(
                            globalOffset: globalOff,
                            batchReads1: subReads1, batchReads2: subReads2,
                            batchRegs1: subRegs1, batchRegs2: subRegs2,
                            doRescue: doRescue, dist: dist, scoring: scoring,
                            gpuEngine: gpuEngine,
                            skipPairing: skipPairing, scoringMat: scoringMat
                        )
                    }
                    activeTasks += 1
                }
            }

            // Drain remaining Phase 3 tasks
            for try await result in phase3Group {
                try handlePhase3Result(result)
            }

            // Print rescue stats (inside task group scope to avoid cross-isolation access)
            let doRescueFinal = insertDist.map {
                !$0.stats[$0.primaryOrientation.rawValue].failed && !skipRescue
            } ?? false
            if doRescueFinal && options.verbosity >= 3 {
                fputs("[PE] Mate rescue: \(rescueCount1) read1 + \(rescueCount2) read2 rescued\n", stderr)
            }
        }

        _ = await phase1Done
    }

    // MARK: - Private Helpers

    /// Promote the paired region to index 0 for primary output.
    /// Matches bwa-mem2's mem_sam_pe() swap logic that ensures the paired region
    /// is emitted as primary. Fixes secondary references accordingly.
    nonisolated private func promotePairedRegion(regions: inout [MemAlnReg], pairedIdx: Int) {
        guard pairedIdx > 0 && pairedIdx < regions.count else { return }

        // bwa-mem2 mem_sam_pe line 456-457: when promoting a secondary, set its sub
        // to the parent primary's score (the region it was originally secondary to).
        // This ensures SE MAPQ reflects the existence of a better-scoring alignment.
        let promoted = regions[pairedIdx]
        if promoted.secondary >= 0 && promoted.secondary < Int32(regions.count) {
            let parentScore = regions[Int(promoted.secondary)].score
            regions[pairedIdx].sub = parentScore
        }

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
        // Use -2 (not -1) to indicate this was a promoted secondary (bwa-mem2 convention)
        regions[0].secondary = -1
    }

    /// Prepare complete output plan for one read pair (parallel-safe).
    /// Handles all PE cases: both-unmapped, pairing resolved, mapped/unmapped, unpaired.
    nonisolated func preparePairOutput(
        read1: ReadSequence, read2: ReadSequence,
        regions1: [MemAlnReg], regions2: [MemAlnReg],
        cigars1: [CIGARInfo], cigars2: [CIGARInfo],
        dist: InsertSizeDistribution,
        genomeLen: Int64,
        skipPairing: Bool,
        scoringMat: [Int8]
    ) -> PairOutputPlan {
        let r1Empty = regions1.isEmpty
        let r2Empty = regions2.isEmpty

        if r1Empty && r2Empty {
            let pe1 = PairedEndInfo(
                isRead1: true, isProperPair: false,
                mateTid: -1, matePos: -1,
                mateIsReverse: false, mateIsUnmapped: true, tlen: 0
            )
            let pe2 = PairedEndInfo(
                isRead1: false, isProperPair: false,
                mateTid: -1, matePos: -1,
                mateIsReverse: false, mateIsUnmapped: true, tlen: 0
            )
            return PairOutputPlan(
                regions1: regions1, regions2: regions2,
                records: [
                    RecordSpec(readIsRead1: true, regionIndex: -1, cigarInfo: nil,
                               mapq: 0, isPrimary: false, isSupplementary: false,
                               isSecondary: false, pairedEnd: pe1, saTag: nil, xaTag: nil),
                    RecordSpec(readIsRead1: false, regionIndex: -1, cigarInfo: nil,
                               mapq: 0, isPrimary: false, isSupplementary: false,
                               isSecondary: false, pairedEnd: pe2, saTag: nil, xaTag: nil),
                ]
            )
        }

        let decision = skipPairing ? nil : PairedEndResolver.resolve(
            regions1: regions1, regions2: regions2,
            dist: dist, genomeLength: genomeLen, scoring: options.scoring
        )

        // bwa-mem2 lines 434-438: check if either read has multiple independent
        // (non-secondary) hits. If so, skip PE pairing and use no_pairing path
        // which outputs supplementary alignments via mem_reg2sam.
        let isMulti1 = Self.hasMultipleNonSecondary(regions1, scoring: options.scoring)
        let isMulti2 = Self.hasMultipleNonSecondary(regions2, scoring: options.scoring)

        if let decision = decision, !isMulti1, !isMulti2 {
            // Both mapped and paired — promote paired regions to primary
            var pairedRegions1 = regions1
            var pairedRegions2 = regions2
            var pairedCigars1 = cigars1
            var pairedCigars2 = cigars2
            promotePairedRegion(regions: &pairedRegions1, pairedIdx: decision.idx1)
            promotePairedRegion(regions: &pairedRegions2, pairedIdx: decision.idx2)

            // Suppress supplementaries: bwa-mem2's PE paired path (lines 489-503)
            // only outputs the paired primary, not supplementary alignments.
            // Mark displaced non-secondary regions as secondary to the promoted primary.
            for i in 1..<pairedRegions1.count {
                if pairedRegions1[i].secondary < 0 {
                    pairedRegions1[i].secondary = 0
                }
            }
            for i in 1..<pairedRegions2.count {
                if pairedRegions2[i].secondary < 0 {
                    pairedRegions2[i].secondary = 0
                }
            }

            if decision.idx1 > 0 && decision.idx1 < pairedCigars1.count {
                pairedCigars1.swapAt(0, decision.idx1)
            }
            if decision.idx2 > 0 && decision.idx2 < pairedCigars2.count {
                pairedCigars2.swapAt(0, decision.idx2)
            }

            // Regenerate CIGAR if a promoted secondary had its CIGAR skipped
            if pairedCigars1[0].cigarString == "*" {
                pairedCigars1[0] = generateCIGAR(read: read1, region: pairedRegions1[0], scoringMatrix: scoringMat)
            }
            if pairedCigars2[0].cigarString == "*" {
                pairedCigars2[0] = generateCIGAR(read: read2, region: pairedRegions2[0], scoringMatrix: scoringMat)
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

            var records: [RecordSpec] = []
            records.append(contentsOf: classifyAlignments(
                read: read1, regions: pairedRegions1, readIsRead1: true,
                pairedEnd: pe1, cigarCache: pairedCigars1, scoringMat: scoringMat
            ))
            records.append(contentsOf: classifyAlignments(
                read: read2, regions: pairedRegions2, readIsRead1: false,
                pairedEnd: pe2, cigarCache: pairedCigars2, scoringMat: scoringMat
            ))

            // Apply PE MAPQ adjustment (bwa-mem2 mem_sam_pe lines 441-467)
            // Find primaries for each read and adjust their MAPQ using pair evidence
            var pri1Idx: Int?
            var pri2Idx: Int?
            for (i, r) in records.enumerated() {
                if r.isPrimary && r.readIsRead1 && pri1Idx == nil { pri1Idx = i }
                if r.isPrimary && !r.readIsRead1 && pri2Idx == nil { pri2Idx = i }
            }
            if let i1 = pri1Idx, let i2 = pri2Idx {
                let (adj1, adj2) = PairedEndResolver.adjustMAPQ(
                    seMAPQ1: records[i1].mapq,
                    seMAPQ2: records[i2].mapq,
                    region1: pairedRegions1[0],
                    region2: pairedRegions2[0],
                    decision: decision,
                    scoring: options.scoring
                )
                records[i1].mapq = adj1
                records[i2].mapq = adj2

                // Cap supplementary MAPQ at adjusted primary MAPQ
                for i in 0..<records.count {
                    if records[i].isSupplementary {
                        if records[i].readIsRead1 {
                            records[i].mapq = min(records[i].mapq, adj1)
                        } else {
                            records[i].mapq = min(records[i].mapq, adj2)
                        }
                    }
                }
            }

            return PairOutputPlan(
                regions1: pairedRegions1, regions2: pairedRegions2, records: records
            )

        } else if !r1Empty && r2Empty {
            return prepareMappedUnmappedPlan(
                mappedRead: read1, mappedRegions: regions1, unmappedRead: read2,
                mappedIsRead1: true, cigarCache: cigars1, scoringMat: scoringMat
            )

        } else if r1Empty && !r2Empty {
            return prepareMappedUnmappedPlan(
                mappedRead: read2, mappedRegions: regions2, unmappedRead: read1,
                mappedIsRead1: false, cigarCache: cigars2, scoringMat: scoringMat
            )

        } else {
            // Both mapped but no valid pairing — use default primaries
            // bwa-mem2 no_pairing (lines 536-541): check if default primaries
            // happen to form a proper pair (same rid, distance within range)
            let cigar1 = cigars1[0]
            let cigar2 = cigars2[0]
            let (rid1, localPos1) = index.metadata.decodePosition(cigar1.pos)
            let (rid2, localPos2) = index.metadata.decodePosition(cigar2.pos)

            var isProper = false
            var tlen1: Int64 = 0
            var tlen2: Int64 = 0

            if rid1 == rid2 && !regions1.isEmpty && !regions2.isEmpty {
                // Check if default primaries form a proper pair
                if let result = InsertSizeEstimator.inferOrientation(
                    r1: regions1[0], r2: regions2[0], genomeLength: genomeLen
                ) {
                    let oriStats = dist.stats[result.orientation.rawValue]
                    if !oriStats.failed
                        && result.insertSize >= oriStats.properLow
                        && result.insertSize <= oriStats.properHigh {
                        isProper = true
                    }
                }
                // Compute TLEN for same-chromosome pairs
                let computed = PairedEndResolver.computeTLEN(
                    pos1: localPos1, isReverse1: cigar1.isReverse, refLen1: cigar1.refConsumed,
                    pos2: localPos2, isReverse2: cigar2.isReverse, refLen2: cigar2.refConsumed
                )
                tlen1 = computed.tlen1
                tlen2 = computed.tlen2
            }

            let pe1 = PairedEndInfo(
                isRead1: true, isProperPair: isProper,
                mateTid: rid2, matePos: localPos2,
                mateIsReverse: cigar2.isReverse, mateIsUnmapped: false,
                tlen: tlen1, mateCigarString: cigar2.cigarString
            )
            let pe2 = PairedEndInfo(
                isRead1: false, isProperPair: isProper,
                mateTid: rid1, matePos: localPos1,
                mateIsReverse: cigar1.isReverse, mateIsUnmapped: false,
                tlen: tlen2, mateCigarString: cigar1.cigarString
            )

            var records: [RecordSpec] = []
            records.append(contentsOf: classifyAlignments(
                read: read1, regions: regions1, readIsRead1: true,
                pairedEnd: pe1, cigarCache: cigars1, scoringMat: scoringMat
            ))
            records.append(contentsOf: classifyAlignments(
                read: read2, regions: regions2, readIsRead1: false,
                pairedEnd: pe2, cigarCache: cigars2, scoringMat: scoringMat
            ))
            return PairOutputPlan(
                regions1: regions1, regions2: regions2, records: records
            )
        }
    }

    /// Prepare output plan for a mapped/unmapped pair (primary only, no supplementaries).
    nonisolated private func prepareMappedUnmappedPlan(
        mappedRead: ReadSequence,
        mappedRegions: [MemAlnReg],
        unmappedRead: ReadSequence,
        mappedIsRead1: Bool,
        cigarCache: [CIGARInfo],
        scoringMat: [Int8]
    ) -> PairOutputPlan {
        let region = mappedRegions[0]
        let cigar = cigarCache.first ?? generateCIGAR(read: mappedRead, region: region, scoringMatrix: scoringMat)
        let (rid, localPos) = index.metadata.decodePosition(cigar.pos)

        let mappedPE = PairedEndInfo(
            isRead1: mappedIsRead1, isProperPair: false,
            mateTid: rid, matePos: localPos,
            mateIsReverse: false, mateIsUnmapped: true, tlen: 0
        )
        let unmappedPE = PairedEndInfo(
            isRead1: !mappedIsRead1, isProperPair: false,
            mateTid: rid, matePos: localPos,
            mateIsReverse: cigar.isReverse, mateIsUnmapped: false, tlen: 0
        )

        let mapq = MappingQuality.compute(
            region: region, allRegions: mappedRegions,
            scoring: options.scoring, readLength: Int32(mappedRead.length)
        )

        var records: [RecordSpec] = []
        let mappedSpec = RecordSpec(
            readIsRead1: mappedIsRead1, regionIndex: 0, cigarInfo: cigar,
            mapq: mapq, isPrimary: true, isSupplementary: false,
            isSecondary: false, pairedEnd: mappedPE, saTag: nil, xaTag: nil
        )
        let unmappedSpec = RecordSpec(
            readIsRead1: !mappedIsRead1, regionIndex: -1, cigarInfo: nil,
            mapq: 0, isPrimary: false, isSupplementary: false,
            isSecondary: false, pairedEnd: unmappedPE, saTag: nil, xaTag: nil
        )

        if mappedIsRead1 {
            records.append(mappedSpec)
            records.append(unmappedSpec)
        } else {
            records.append(unmappedSpec)
            records.append(mappedSpec)
        }

        let (r1, r2) = mappedIsRead1
            ? (mappedRegions, [MemAlnReg]())
            : ([MemAlnReg](), mappedRegions)
        return PairOutputPlan(regions1: r1, regions2: r2, records: records)
    }

    /// Process a batch of pairs: rescue + CIGAR + output preparation (fully merged).
    /// When gpuEngine is available, collects all rescue candidates for the batch,
    /// dispatches GPU batches, applies results, generates CIGARs, then builds output plans.
    nonisolated private func processFullBatch(
        globalOffset: Int,
        batchReads1: [ReadSequence], batchReads2: [ReadSequence],
        batchRegs1: [[MemAlnReg]], batchRegs2: [[MemAlnReg]],
        doRescue: Bool, dist: InsertSizeDistribution,
        scoring: ScoringParameters,
        gpuEngine: Any?,
        skipPairing: Bool, scoringMat: [Int8]
    ) async -> [(Int, PairOutputPlan, Int, Int)] {
        let mat = scoring.scoringMatrix()
        let count = batchReads1.count
        var finalRegs1 = batchRegs1
        var finalRegs2 = batchRegs2
        var rescueCounts1 = [Int](repeating: 0, count: count)
        var rescueCounts2 = [Int](repeating: 0, count: count)

        if doRescue {
            #if canImport(Metal)
            if let engine = gpuEngine as? MetalSWEngine {
                // GPU per-batch rescue: collect candidates, dispatch, apply
                let genomeLen = index.genomeLength

                // --- Collect ALL rescue candidates (both directions) upfront ---
                var allCands: [MateRescue.RescueCandidate] = []
                var candRanges2 = [(start: Int, count: Int)](repeating: (0, 0), count: count)
                var candRanges1 = [(start: Int, count: Int)](repeating: (0, 0), count: count)

                // Direction 1: rescue R2 from R1 templates
                for li in 0..<count {
                    let templates = Self.selectRescueCandidates(finalRegs1[li], scoring: scoring)
                    guard !templates.isEmpty else { continue }
                    let cands = MateRescue.collectCandidates(
                        templateRegions: templates,
                        mateRead: batchReads2[li],
                        mateRegions: finalRegs2[li],
                        dist: dist,
                        genomeLength: genomeLen,
                        packedRef: index.packedRef,
                        metadata: index.metadata
                    )
                    if !cands.isEmpty {
                        candRanges2[li] = (start: allCands.count, count: cands.count)
                        allCands.append(contentsOf: cands)
                    }
                }

                // Direction 2: rescue R1 from original R2 templates
                for li in 0..<count {
                    let templates = Self.selectRescueCandidates(finalRegs2[li], scoring: scoring)
                    guard !templates.isEmpty else { continue }
                    let cands = MateRescue.collectCandidates(
                        templateRegions: templates,
                        mateRead: batchReads1[li],
                        mateRegions: finalRegs1[li],
                        dist: dist,
                        genomeLength: genomeLen,
                        packedRef: index.packedRef,
                        metadata: index.metadata
                    )
                    if !cands.isEmpty {
                        candRanges1[li] = (start: allCands.count, count: cands.count)
                        allCands.append(contentsOf: cands)
                    }
                }

                // --- Single GPU dispatch for both directions ---
                if !allCands.isEmpty {
                    let tasks = allCands.map {
                        LocalSWTask(query: $0.mateSeq, target: $0.refBases, scoring: scoring)
                    }
                    let gpuResults = await LocalSWDispatcher.dispatchBatch(tasks: tasks, engine: engine)

                    let minSubScore = UInt8(clamping: scoring.minSeedLength &* scoring.matchScore)

                    // Apply direction 1 results: rescue R2
                    for li in 0..<count {
                        let r = candRanges2[li]
                        for j in 0..<r.count {
                            let ci = r.start + j
                            let cand = allCands[ci]
                            let sw: (score: Int32, qb: Int32, qe: Int32, tb: Int32, te: Int32, sc2: Int32)?
                            if let g = gpuResults[ci] {
                                // GPU path doesn't compute score2 yet
                                sw = (g.score, g.queryBegin, g.queryEnd, g.targetBegin, g.targetEnd, 0)
                            } else if let c = LocalSWAligner.align(
                                query: cand.mateSeq, target: cand.refBases,
                                scoring: scoring, scoringMatrix: mat,
                                minSubScore: minSubScore
                            ) {
                                sw = (c.score, c.queryBegin, c.queryEnd, c.targetBegin, c.targetEnd, c.score2)
                            } else { sw = nil }
                            if let sw = sw, let reg = MateRescue.applyResult(
                                candidate: cand, score: sw.score,
                                queryBegin: sw.qb, queryEnd: sw.qe,
                                targetBegin: sw.tb, targetEnd: sw.te,
                                score2: sw.sc2,
                                genomeLength: genomeLen, scoring: scoring
                            ) {
                                finalRegs2[li].append(reg)
                                rescueCounts2[li] += 1
                            }
                        }
                    }

                    // Apply direction 2 results: rescue R1
                    for li in 0..<count {
                        let r = candRanges1[li]
                        for j in 0..<r.count {
                            let ci = r.start + j
                            let cand = allCands[ci]
                            let sw: (score: Int32, qb: Int32, qe: Int32, tb: Int32, te: Int32, sc2: Int32)?
                            if let g = gpuResults[ci] {
                                sw = (g.score, g.queryBegin, g.queryEnd, g.targetBegin, g.targetEnd, 0)
                            } else if let c = LocalSWAligner.align(
                                query: cand.mateSeq, target: cand.refBases,
                                scoring: scoring, scoringMatrix: mat,
                                minSubScore: minSubScore
                            ) {
                                sw = (c.score, c.queryBegin, c.queryEnd, c.targetBegin, c.targetEnd, c.score2)
                            } else { sw = nil }
                            if let sw = sw, let reg = MateRescue.applyResult(
                                candidate: cand, score: sw.score,
                                queryBegin: sw.qb, queryEnd: sw.qe,
                                targetBegin: sw.tb, targetEnd: sw.te,
                                score2: sw.sc2,
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
                    let gi = globalOffset + li
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
                    let gi = globalOffset + li
                    let (_, rr1, rr2, c1, c2) = rescuePair(
                        idx: gi, regions1: finalRegs1[li], regions2: finalRegs2[li],
                        read1: batchReads1[li], read2: batchReads2[li], dist: dist, scoring: scoring
                    )
                    finalRegs1[li] = rr1; finalRegs2[li] = rr2
                    rescueCounts1[li] = c1; rescueCounts2[li] = c2
                }
            }
            #else
            for li in 0..<count {
                let gi = globalOffset + li
                let (_, rr1, rr2, c1, c2) = rescuePair(
                    idx: gi, regions1: finalRegs1[li], regions2: finalRegs2[li],
                    read1: batchReads1[li], read2: batchReads2[li], dist: dist, scoring: scoring
                )
                finalRegs1[li] = rr1; finalRegs2[li] = rr2
                rescueCounts1[li] = c1; rescueCounts2[li] = c2
            }
            #endif
        }

        // CIGAR generation + output preparation
        var results: [(Int, PairOutputPlan, Int, Int)] = []
        results.reserveCapacity(count)
        for li in 0..<count {
            let gi = globalOffset + li
            let cig1 = generateFilteredCIGARs(read: batchReads1[li], regions: finalRegs1[li], scoringMatrix: mat)
            let cig2 = generateFilteredCIGARs(read: batchReads2[li], regions: finalRegs2[li], scoringMatrix: mat)
            let plan = preparePairOutput(
                read1: batchReads1[li], read2: batchReads2[li],
                regions1: finalRegs1[li], regions2: finalRegs2[li],
                cigars1: cig1, cigars2: cig2,
                dist: dist, genomeLen: index.genomeLength,
                skipPairing: skipPairing, scoringMat: scoringMat
            )
            results.append((gi, plan, rescueCounts1[li], rescueCounts2[li]))
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

    /// Check if regions have multiple independent (non-secondary) hits above score threshold.
    /// Matches bwa-mem2's is_multi check (bwamem_pair.cpp lines 434-438).
    private static func hasMultipleNonSecondary(
        _ regions: [MemAlnReg], scoring: ScoringParameters
    ) -> Bool {
        var count = 0
        for r in regions where r.secondary < 0 && r.score >= scoring.minOutputScore {
            count += 1
            if count > 1 { return true }
        }
        return false
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

        // Double-buffered loop: overlap GPU SMEM of batch N+1 with CPU extension of batch N.
        // Prime the pipeline by dispatching the first batch to GPU (non-blocking).
        var batchStart = 0
        var batchEnd = min(batchStart + gpuBatchSize, readCount)
        var batchReads = Array(reads[batchStart..<batchEnd])
        var queries = batchReads.map { $0.bases }
        var handle = SMEMDispatcher.dispatchBatchAsync(queries: queries, engine: engine)

        while batchStart < readCount {
            let batchCount = batchEnd - batchStart

            // Wait for GPU results of current batch
            let fwdResults: [ForwardExtResult]
            if let h = handle {
                fwdResults = SMEMDispatcher.collectResults(handle: h)
            } else {
                fwdResults = queries.map { q in
                    ForwardExtResult(
                        rightEnds: [Int32](repeating: 0, count: q.count),
                        kValues: [Int64](repeating: 0, count: q.count),
                        sValues: [Int64](repeating: 0, count: q.count)
                    )
                }
            }
            let gpuSMEMs = SMEMDispatcher.extractSMEMs(
                from: fwdResults,
                queries: queries,
                minSeedLen: scoring.minSeedLength
            )

            // Immediately dispatch NEXT batch to GPU (non-blocking) so it runs
            // concurrently with the CPU extension of the current batch below.
            let nextBatchStart = batchEnd
            let nextBatchEnd = min(nextBatchStart + gpuBatchSize, readCount)
            var nextHandle: SMEMDispatchHandle? = nil
            var nextBatchReads: [ReadSequence] = []
            var nextQueries: [[UInt8]] = []
            if nextBatchStart < readCount {
                nextBatchReads = Array(reads[nextBatchStart..<nextBatchEnd])
                nextQueries = nextBatchReads.map { $0.bases }
                nextHandle = SMEMDispatcher.dispatchBatchAsync(
                    queries: nextQueries, engine: engine)
            }

            // Phases 1.5-5: parallel per-read on CPU (overlaps with GPU running next batch)
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

            // Advance to next batch
            batchStart = nextBatchStart
            batchEnd = nextBatchEnd
            batchReads = nextBatchReads
            queries = nextQueries
            handle = nextHandle
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

    /// CPU streaming alignment: processes reads in 4096-read sub-batches,
    /// yielding CompletedBatch results through an AsyncStream continuation
    /// so downstream Phase 3 can begin before all reads are aligned.
    nonisolated private func alignReadsCPUBatched(
        _ reads: [ReadSequence],
        mate: Int,
        idBase: UInt64, idShift: Bool,
        continuation: AsyncStream<CompletedBatch>.Continuation
    ) async {
        let maxConcurrency = options.scoring.numThreads
        let cpuBatchSize = 4096
        let readCount = reads.count
        var batchIdx = 0

        for batchStart in stride(from: 0, to: readCount, by: cpuBatchSize) {
            let batchEnd = min(batchStart + cpuBatchSize, readCount)
            let batchCount = batchEnd - batchStart
            let batchReads = Array(reads[batchStart..<batchEnd])

            let results = await withTaskGroup(
                of: (Int, [MemAlnReg]).self,
                returning: [(Int, [MemAlnReg])].self
            ) { group in
                var nextIdx = 0
                var collected: [(Int, [MemAlnReg])] = []
                collected.reserveCapacity(batchCount)

                while nextIdx < min(maxConcurrency, batchCount) {
                    let idx = nextIdx
                    let read = batchReads[idx]
                    let gi = batchStart + idx
                    let readId = idShift ? (UInt64(gi) &<< 1) | idBase : UInt64(gi)
                    group.addTask { [self] in
                        (idx, self.alignRead(read, readId: readId))
                    }
                    nextIdx += 1
                }

                for await result in group {
                    collected.append(result)
                    if nextIdx < batchCount {
                        let idx = nextIdx
                        let read = batchReads[idx]
                        let gi = batchStart + idx
                        let readId = idShift ? (UInt64(gi) &<< 1) | idBase : UInt64(gi)
                        group.addTask { [self] in
                            (idx, self.alignRead(read, readId: readId))
                        }
                        nextIdx += 1
                    }
                }

                return collected
            }

            var regions = Array(repeating: [MemAlnReg](), count: batchCount)
            for (idx, regs) in results {
                regions[idx] = regs
            }
            continuation.yield(CompletedBatch(
                mate: mate, batchIndex: batchIdx,
                batchStart: batchStart,
                reads: batchReads, regions: regions
            ))
            batchIdx += 1
        }
    }

    #if canImport(Metal)
    /// GPU-seeded streaming alignment: same double-buffered GPU SMEM loop as
    /// alignAllReadsGPUSeeded, but yields CompletedBatch results through an
    /// AsyncStream continuation instead of accumulating into a local array.
    nonisolated private func alignReadsGPUSeededStreaming(
        _ reads: [ReadSequence],
        mate: Int,
        idBase: UInt64, idShift: Bool,
        engine: MetalSWEngine,
        continuation: AsyncStream<CompletedBatch>.Continuation
    ) async {
        let scoring = options.scoring
        let maxConcurrency = scoring.numThreads
        let gpuBatchSize = 4096
        let readCount = reads.count

        // Double-buffered loop: overlap GPU SMEM of batch N+1 with CPU extension of batch N.
        var batchStart = 0
        var batchEnd = min(batchStart + gpuBatchSize, readCount)
        var batchReads = Array(reads[batchStart..<batchEnd])
        var queries = batchReads.map { $0.bases }
        var handle = SMEMDispatcher.dispatchBatchAsync(queries: queries, engine: engine)
        var batchIdx = 0

        while batchStart < readCount {
            let batchCount = batchEnd - batchStart

            // Wait for GPU results of current batch
            let fwdResults: [ForwardExtResult]
            if let h = handle {
                fwdResults = SMEMDispatcher.collectResults(handle: h)
            } else {
                fwdResults = queries.map { q in
                    ForwardExtResult(
                        rightEnds: [Int32](repeating: 0, count: q.count),
                        kValues: [Int64](repeating: 0, count: q.count),
                        sValues: [Int64](repeating: 0, count: q.count)
                    )
                }
            }
            let gpuSMEMs = SMEMDispatcher.extractSMEMs(
                from: fwdResults,
                queries: queries,
                minSeedLen: scoring.minSeedLength
            )

            // Immediately dispatch NEXT batch to GPU (non-blocking)
            let nextBatchStart = batchEnd
            let nextBatchEnd = min(nextBatchStart + gpuBatchSize, readCount)
            var nextHandle: SMEMDispatchHandle? = nil
            var nextBatchReads: [ReadSequence] = []
            var nextQueries: [[UInt8]] = []
            if nextBatchStart < readCount {
                nextBatchReads = Array(reads[nextBatchStart..<nextBatchEnd])
                nextQueries = nextBatchReads.map { $0.bases }
                nextHandle = SMEMDispatcher.dispatchBatchAsync(
                    queries: nextQueries, engine: engine)
            }

            // Phases 1.5-5: parallel per-read on CPU (overlaps with GPU running next batch)
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

            // Yield completed sub-batch
            var regions = Array(repeating: [MemAlnReg](), count: batchCount)
            for (idx, regs) in batchResults {
                regions[idx] = regs
            }
            continuation.yield(CompletedBatch(
                mate: mate, batchIndex: batchIdx,
                batchStart: capturedBatchStart,
                reads: batchReads, regions: regions
            ))

            // Advance to next batch
            batchStart = nextBatchStart
            batchEnd = nextBatchEnd
            batchReads = nextBatchReads
            queries = nextQueries
            handle = nextHandle
            batchIdx += 1
        }
    }
    #endif

    /// Launch Phase 1 streaming producers for both mates concurrently,
    /// then finish the continuation when both complete.
    nonisolated private func runPhase1Streaming(
        reads1: [ReadSequence], reads2: [ReadSequence],
        gpuEngine: Any?,
        continuation: AsyncStream<CompletedBatch>.Continuation
    ) async {
        await withTaskGroup(of: Void.self) { group in
            #if canImport(Metal)
            if let engine = gpuEngine as? MetalSWEngine {
                group.addTask { [self] in
                    await self.alignReadsGPUSeededStreaming(
                        reads1, mate: 1, idBase: 0, idShift: true,
                        engine: engine, continuation: continuation
                    )
                }
                group.addTask { [self] in
                    await self.alignReadsGPUSeededStreaming(
                        reads2, mate: 2, idBase: 1, idShift: true,
                        engine: engine, continuation: continuation
                    )
                }
            } else {
                group.addTask { [self] in
                    await self.alignReadsCPUBatched(
                        reads1, mate: 1, idBase: 0, idShift: true,
                        continuation: continuation
                    )
                }
                group.addTask { [self] in
                    await self.alignReadsCPUBatched(
                        reads2, mate: 2, idBase: 1, idShift: true,
                        continuation: continuation
                    )
                }
            }
            #else
            group.addTask { [self] in
                await self.alignReadsCPUBatched(
                    reads1, mate: 1, idBase: 0, idShift: true,
                    continuation: continuation
                )
            }
            group.addTask { [self] in
                await self.alignReadsCPUBatched(
                    reads2, mate: 2, idBase: 1, idShift: true,
                    continuation: continuation
                )
            }
            #endif
        }
        continuation.finish()
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

}
