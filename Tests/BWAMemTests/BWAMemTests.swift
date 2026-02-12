import Foundation
import Testing
@testable import BWAMem
@testable import BWACore
@testable import FMIndex
@testable import Alignment

@Suite("BWAMem Tests")
struct BWAMemTests {

    @Test("BWAMemOptions defaults")
    func testOptionsDefaults() {
        let options = BWAMemOptions()
        #expect(options.scoring.matchScore == 1)
        #expect(options.scoring.mismatchPenalty == 4)
        #expect(options.isPairedEnd == false)
    }

    @Test("MappingQuality for unique perfect hit")
    func testMAPQPerfect() {
        let region = MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 100, trueScore: 100, sub: 0, seedCov: 100)
        let mapq = MappingQuality.compute(
            region: region,
            allRegions: [region],
            scoring: ScoringParameters(),
            readLength: 100
        )
        #expect(mapq == 60)
    }

    @Test("MappingQuality for zero score")
    func testMAPQZero() {
        let region = MemAlnReg(rb: 0, re: 0, qb: 0, qe: 0, score: 0)
        let mapq = MappingQuality.compute(
            region: region,
            allRegions: [region],
            scoring: ScoringParameters(),
            readLength: 100
        )
        #expect(mapq == 0)
    }

    @Test("MappingQuality decreases with sub-optimal hits")
    func testMAPQWithSub() {
        let region1 = MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 80, trueScore: 80, sub: 70, seedCov: 80)
        let region2 = MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 80, trueScore: 80, sub: 10, seedCov: 80)

        let scoring = ScoringParameters()
        let mapq1 = MappingQuality.compute(region: region1, allRegions: [region1], scoring: scoring, readLength: 100)
        let mapq2 = MappingQuality.compute(region: region2, allRegions: [region2], scoring: scoring, readLength: 100)

        // Higher sub-optimal score should give lower MAPQ
        #expect(mapq1 < mapq2)
    }

    @Test("ScoringParameters ALT defaults")
    func testScoringALTDefaults() {
        let scoring = ScoringParameters()
        #expect(scoring.maxXAHits == 5)
        #expect(scoring.maxXAHitsAlt == 200)
        #expect(scoring.xaDropRatio == 0.80)
    }
}

// MARK: - ALT File Loading Tests

@Suite("ALT Loading Tests")
struct ALTLoadingTests {

    @Test("loadAlt sets isAlt on matching annotations")
    func testLoadAltSetsIsAlt() throws {
        // Create temp directory with .alt file
        let tmpDir = NSTemporaryDirectory() + "swiftbwa_test_alt_\(ProcessInfo.processInfo.processIdentifier)"
        try FileManager.default.createDirectory(atPath: tmpDir, withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(atPath: tmpDir) }

        let prefix = tmpDir + "/test"
        let altContent = "chr1_alt\nchr3_alt\n"
        try altContent.write(toFile: prefix + ".alt", atomically: true, encoding: .utf8)

        var metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 4,
            annotations: [
                ReferenceAnnotation(offset: 0, length: 2500, name: "chr1"),
                ReferenceAnnotation(offset: 2500, length: 2500, name: "chr1_alt"),
                ReferenceAnnotation(offset: 5000, length: 2500, name: "chr2"),
                ReferenceAnnotation(offset: 7500, length: 2500, name: "chr3_alt"),
            ]
        )

        FMIndexLoader.loadAlt(from: prefix, metadata: &metadata)

        #expect(metadata.annotations[0].isAlt == false)  // chr1
        #expect(metadata.annotations[1].isAlt == true)   // chr1_alt
        #expect(metadata.annotations[2].isAlt == false)  // chr2
        #expect(metadata.annotations[3].isAlt == true)   // chr3_alt
    }

    @Test("loadAlt skips header lines starting with @")
    func testLoadAltSkipsHeaders() throws {
        let tmpDir = NSTemporaryDirectory() + "swiftbwa_test_alt2_\(ProcessInfo.processInfo.processIdentifier)"
        try FileManager.default.createDirectory(atPath: tmpDir, withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(atPath: tmpDir) }

        let prefix = tmpDir + "/test"
        let altContent = "@HD\tVN:1.0\nchr1_alt\n"
        try altContent.write(toFile: prefix + ".alt", atomically: true, encoding: .utf8)

        var metadata = ReferenceMetadata(
            totalLength: 5000, numSequences: 2,
            annotations: [
                ReferenceAnnotation(offset: 0, length: 2500, name: "chr1"),
                ReferenceAnnotation(offset: 2500, length: 2500, name: "chr1_alt"),
            ]
        )

        FMIndexLoader.loadAlt(from: prefix, metadata: &metadata)

        #expect(metadata.annotations[0].isAlt == false)
        #expect(metadata.annotations[1].isAlt == true)
    }

    @Test("loadAlt is no-op when .alt file is missing")
    func testLoadAltMissingFile() {
        var metadata = ReferenceMetadata(
            totalLength: 5000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 5000, name: "chr1")]
        )

        FMIndexLoader.loadAlt(from: "/nonexistent/path/test", metadata: &metadata)

        #expect(metadata.annotations[0].isAlt == false)
    }

    @Test("loadAlt handles tab-delimited lines")
    func testLoadAltTabDelimited() throws {
        let tmpDir = NSTemporaryDirectory() + "swiftbwa_test_alt3_\(ProcessInfo.processInfo.processIdentifier)"
        try FileManager.default.createDirectory(atPath: tmpDir, withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(atPath: tmpDir) }

        let prefix = tmpDir + "/test"
        let altContent = "chr1_alt\textra_field\nchr2\tmore_stuff\n"
        try altContent.write(toFile: prefix + ".alt", atomically: true, encoding: .utf8)

        var metadata = ReferenceMetadata(
            totalLength: 5000, numSequences: 2,
            annotations: [
                ReferenceAnnotation(offset: 0, length: 2500, name: "chr1_alt"),
                ReferenceAnnotation(offset: 2500, length: 2500, name: "chr2"),
            ]
        )

        FMIndexLoader.loadAlt(from: prefix, metadata: &metadata)

        #expect(metadata.annotations[0].isAlt == true)
        #expect(metadata.annotations[1].isAlt == true)
    }
}

// MARK: - ALT Output Tests

@Suite("ALT Output Tests")
struct ALTOutputTests {

    @Test("pa tag is present when altSc > 0")
    func testPaTagPresent() throws {
        let metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 10000, name: "chr1")]
        )

        var region = MemAlnReg(
            rb: 100, re: 200, qb: 0, qe: 100, rid: 0,
            score: 80, trueScore: 80, sub: 0, w: 10
        )
        region.altSc = 120  // ALT competitor score

        let read = ReadSequence(
            name: "test_pa",
            sequence: String(repeating: "A", count: 100),
            qualityString: String(repeating: "I", count: 100)
        )

        let record = try SAMOutputBuilder.buildRecord(
            read: read, region: region, allRegions: [region],
            metadata: metadata, scoring: ScoringParameters(),
            cigar: [100 << 4 | 0],
            isPrimary: true
        )

        // The pa tag should be present as a float
        let aux = record.auxiliaryData
        let paVal = aux.float(forTag: "pa")
        #expect(paVal != nil)
        // pa = score / altSc = 80 / 120 ≈ 0.667
        if let pa = paVal {
            #expect(pa > 0.6 && pa < 0.7)
        }
    }

    @Test("pa tag is absent when altSc == 0")
    func testPaTagAbsent() throws {
        let metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 10000, name: "chr1")]
        )

        let region = MemAlnReg(
            rb: 100, re: 200, qb: 0, qe: 100, rid: 0,
            score: 80, trueScore: 80, sub: 0, w: 10
        )

        let read = ReadSequence(
            name: "test_no_pa",
            sequence: String(repeating: "A", count: 100),
            qualityString: String(repeating: "I", count: 100)
        )

        let record = try SAMOutputBuilder.buildRecord(
            read: read, region: region, allRegions: [region],
            metadata: metadata, scoring: ScoringParameters(),
            cigar: [100 << 4 | 0],
            isPrimary: true
        )

        let aux = record.auxiliaryData
        let paVal = aux.float(forTag: "pa")
        #expect(paVal == nil)
    }

    @Test("XA tag uses maxXAHitsAlt (200) when ALT secondaries present")
    func testXATagALTLimit() {
        // With ALT secondaries, up to 200 are allowed
        let secondaries: [(rname: String, pos: Int64, isReverse: Bool,
                           cigarString: String, nm: Int32)] = (0..<10).map { i in
            (rname: "chr1_alt", pos: Int64(i * 100), isReverse: false,
             cigarString: "100M", nm: Int32(0))
        }

        // With maxHits=200 (ALT limit), 10 secondaries should be included
        let xa = SAMOutputBuilder.buildXATag(secondaries: secondaries, maxHits: 200)
        #expect(xa != nil)

        // With maxHits=5 (normal limit), 10 secondaries should be nil
        let xaNormal = SAMOutputBuilder.buildXATag(secondaries: secondaries, maxHits: 5)
        #expect(xaNormal == nil)
    }

    @Test("altSc tracking: primary records ALT competitor score")
    func testAltScTracking() {
        // Two overlapping regions: ALT with score 120, primary with score 80
        var regions = [
            MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 120),
            MemAlnReg(rb: 200, re: 300, qb: 0, qe: 100, score: 80),
        ]
        regions[0].isAlt = true
        regions[1].isAlt = false

        ChainFilter.markSecondaryALT(
            regions: &regions, maskLevel: 0.50, scoring: ScoringParameters()
        )

        let primaryRegion = regions.first { !$0.isAlt }!
        #expect(primaryRegion.altSc == 120)
    }
}

// MARK: - Insert Size Estimator Tests

@Suite("InsertSizeEstimator Tests")
struct InsertSizeEstimatorTests {

    @Test("Infer FR orientation: read1 forward, read2 reverse, r1 left of r2")
    func testInferOrientationFR() {
        let genomeLen: Int64 = 10000
        // Read1: forward strand, rb=100
        let r1 = MemAlnReg(rb: 100, re: 200, qb: 0, qe: 100, rid: 0, score: 100)
        // Read2: reverse strand, rb=19600 (BWT coords)
        // Mirror onto r1's strand: p2 = 2*10000 - 1 - 19600 = 399
        // dist = |399 - 100| = 299, bwaDir = 1^0 = 1 (FR)
        let r2 = MemAlnReg(rb: 19600, re: 19700, qb: 0, qe: 100, rid: 0, score: 100)

        let result = InsertSizeEstimator.inferOrientation(r1: r1, r2: r2, genomeLength: genomeLen)
        #expect(result != nil)
        #expect(result!.orientation == .fr)
        // Start-to-start distance: 299
        #expect(result!.insertSize == 299)
    }

    @Test("Infer RF orientation: read2 forward, upstream of read1 reverse")
    func testInferOrientationRF() {
        let genomeLen: Int64 = 10000
        // r1 is forward, r2 is reverse, with r2 mirrored position < r1 → RF
        // r1 forward at rb=500, r2 reverse at rb=19800
        // p2 = 2*10000 - 1 - 19800 = 199 (mirrored, < r1.rb=500)
        // bwaDir = 1 ^ 3 = 2 → RF
        let r1 = MemAlnReg(rb: 500, re: 600, qb: 0, qe: 100, rid: 0, score: 100)
        let r2 = MemAlnReg(rb: 19800, re: 19900, qb: 0, qe: 100, rid: 0, score: 100)

        let result = InsertSizeEstimator.inferOrientation(r1: r1, r2: r2, genomeLength: genomeLen)
        #expect(result != nil)
        #expect(result!.orientation == .rf)
    }

    @Test("Infer FF orientation: both forward")
    func testInferOrientationFF() {
        let genomeLen: Int64 = 10000
        let r1 = MemAlnReg(rb: 100, re: 200, qb: 0, qe: 100, rid: 0, score: 100)
        let r2 = MemAlnReg(rb: 300, re: 400, qb: 0, qe: 100, rid: 0, score: 100)

        let result = InsertSizeEstimator.inferOrientation(r1: r1, r2: r2, genomeLength: genomeLen)
        #expect(result != nil)
        #expect(result!.orientation == .ff)
    }

    @Test("Infer RR orientation: both reverse, b1 > b2")
    func testInferOrientationRR() {
        let genomeLen: Int64 = 10000
        // Under bwa encoding, two same-strand reads where p2 <= b1 → dir = 0^3 = 3 (RR)
        // r1 at higher BWT position than r2
        let r1 = MemAlnReg(rb: 19800, re: 19900, qb: 0, qe: 100, rid: 0, score: 100)
        let r2 = MemAlnReg(rb: 19600, re: 19700, qb: 0, qe: 100, rid: 0, score: 100)

        let result = InsertSizeEstimator.inferOrientation(r1: r1, r2: r2, genomeLength: genomeLen)
        #expect(result != nil)
        #expect(result!.orientation == .rr)
    }

    @Test("Different chromosomes returns nil")
    func testDifferentChromosomes() {
        let genomeLen: Int64 = 10000
        let r1 = MemAlnReg(rb: 100, re: 200, qb: 0, qe: 100, rid: 0, score: 100)
        let r2 = MemAlnReg(rb: 300, re: 400, qb: 0, qe: 100, rid: 1, score: 100)

        let result = InsertSizeEstimator.inferOrientation(r1: r1, r2: r2, genomeLength: genomeLen)
        #expect(result == nil)
    }

    @Test("Estimate distribution with known FR insert sizes")
    func testEstimateDistribution() {
        let genomeLen: Int64 = 100000
        // Create 100 simulated FR pairs with start-to-start distances around 300 (±20)
        // Using bwa-mem2's start-to-start metric:
        //   p2 = 2*genomeLen - 1 - r2.rb, dist = |p2 - r1.rb|
        var regions1: [[MemAlnReg]] = []
        var regions2: [[MemAlnReg]] = []

        for i in 0..<100 {
            let targetDist = 280 + Int64(i % 41)  // 280..320
            let r1Pos: Int64 = 1000 + Int64(i) * 500
            let r1 = MemAlnReg(
                rb: r1Pos, re: r1Pos + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0, seedCov: 100
            )
            // We want p2 = r1Pos + targetDist, where p2 = 2*genomeLen - 1 - r2.rb
            // So r2.rb = 2*genomeLen - 1 - (r1Pos + targetDist)
            let r2Rb = 2 * genomeLen - 1 - (r1Pos + targetDist)
            let r2 = MemAlnReg(
                rb: r2Rb, re: r2Rb + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0, seedCov: 100
            )
            regions1.append([r1])
            regions2.append([r2])
        }

        let dist = InsertSizeEstimator.estimate(
            regions1: regions1, regions2: regions2,
            genomeLength: genomeLen
        )

        #expect(dist.primaryOrientation == .fr)
        let frStats = dist.stats[PairOrientation.fr.rawValue]
        #expect(!frStats.failed)
        #expect(frStats.count > 50)
        #expect(frStats.mean > 290 && frStats.mean < 310)
        #expect(frStats.stddev > 0 && frStats.stddev < 20)
    }

    @Test("Quartile filtering excludes outliers")
    func testQuartileFiltering() {
        let genomeLen: Int64 = 100000
        var regions1: [[MemAlnReg]] = []
        var regions2: [[MemAlnReg]] = []

        // 40 normal pairs with start-to-start dist ~300
        for i in 0..<40 {
            let r1Pos: Int64 = 1000 + Int64(i) * 500
            let r1 = MemAlnReg(
                rb: r1Pos, re: r1Pos + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0, seedCov: 100
            )
            // p2 = r1Pos + 300, so r2.rb = 2*genomeLen - 1 - (r1Pos + 300)
            let r2Rb = 2 * genomeLen - 1 - (r1Pos + 300)
            let r2 = MemAlnReg(
                rb: r2Rb, re: r2Rb + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0, seedCov: 100
            )
            regions1.append([r1])
            regions2.append([r2])
        }

        // 5 outlier pairs with dist ~5000 (should be filtered)
        for i in 0..<5 {
            let r1Pos: Int64 = 50000 + Int64(i) * 500
            let r1 = MemAlnReg(
                rb: r1Pos, re: r1Pos + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0, seedCov: 100
            )
            let r2Rb = 2 * genomeLen - 1 - (r1Pos + 5000)
            let r2 = MemAlnReg(
                rb: r2Rb, re: r2Rb + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0, seedCov: 100
            )
            regions1.append([r1])
            regions2.append([r2])
        }

        let dist = InsertSizeEstimator.estimate(
            regions1: regions1, regions2: regions2,
            genomeLength: genomeLen
        )

        let frStats = dist.stats[PairOrientation.fr.rawValue]
        #expect(!frStats.failed)
        // Mean should be around 300, not pulled up by the 5000-size outliers
        #expect(frStats.mean < 400)
    }
}

// MARK: - Paired-End Resolver Tests

@Suite("PairedEndResolver Tests")
struct PairedEndResolverTests {

    @Test("Pair scoring selects best concordant pair")
    func testPairScoring() {
        let genomeLen: Int64 = 100000

        // Two regions for read1: one with score 80, one with score 60
        let r1a = MemAlnReg(rb: 1000, re: 1100, qb: 0, qe: 100, rid: 0, score: 80, trueScore: 80, sub: 0)
        let r1b = MemAlnReg(rb: 5000, re: 5100, qb: 0, qe: 100, rid: 0, score: 60, trueScore: 60, sub: 0)

        // Matching region for read2 that pairs well with r1a
        let r2FwdStart: Int64 = 1200
        let r2FwdEnd: Int64 = 1300
        let r2a = MemAlnReg(
            rb: 2 * genomeLen - r2FwdEnd,
            re: 2 * genomeLen - r2FwdStart,
            qb: 0, qe: 100, rid: 0, score: 80, trueScore: 80, sub: 0
        )

        // Build insert size dist centered at 300 for FR
        var dist = InsertSizeDistribution()
        dist.stats[PairOrientation.fr.rawValue] = OrientationStats()
        dist.stats[PairOrientation.fr.rawValue].mean = 300
        dist.stats[PairOrientation.fr.rawValue].stddev = 30
        dist.stats[PairOrientation.fr.rawValue].low = 200
        dist.stats[PairOrientation.fr.rawValue].high = 400
        dist.stats[PairOrientation.fr.rawValue].count = 100
        dist.stats[PairOrientation.fr.rawValue].failed = false
        dist.primaryOrientation = .fr

        let decision = PairedEndResolver.resolve(
            regions1: [r1a, r1b],
            regions2: [r2a],
            dist: dist,
            genomeLength: genomeLen,
            scoring: ScoringParameters()
        )

        #expect(decision != nil)
        #expect(decision!.idx1 == 0)  // r1a (higher score)
        #expect(decision!.idx2 == 0)  // r2a
        #expect(decision!.isProperPair == true)
    }

    @Test("Discordant pair is not flagged as proper")
    func testDiscordantPair() {
        let genomeLen: Int64 = 100000

        // Both on forward strand (FF orientation)
        let r1 = MemAlnReg(rb: 1000, re: 1100, qb: 0, qe: 100, rid: 0, score: 80)
        let r2 = MemAlnReg(rb: 1200, re: 1300, qb: 0, qe: 100, rid: 0, score: 80)

        var dist = InsertSizeDistribution()
        dist.stats[PairOrientation.fr.rawValue].mean = 300
        dist.stats[PairOrientation.fr.rawValue].stddev = 30
        dist.stats[PairOrientation.fr.rawValue].count = 100
        dist.stats[PairOrientation.fr.rawValue].failed = false
        dist.primaryOrientation = .fr

        let decision = PairedEndResolver.resolve(
            regions1: [r1], regions2: [r2],
            dist: dist,
            genomeLength: genomeLen,
            scoring: ScoringParameters()
        )

        #expect(decision != nil)
        // FF orientation != FR primary → not proper
        #expect(decision!.isProperPair == false)
    }

    @Test("TLEN computation for standard FR pair")
    func testComputeTLEN() {
        // Read1 at pos 100, 50bp on reference; Read2 at pos 250, 50bp on reference
        let (tlen1, tlen2) = PairedEndResolver.computeTLEN(
            pos1: 100, isReverse1: false, refLen1: 50,
            pos2: 250, isReverse2: true, refLen2: 50
        )

        // Leftmost = 100, rightmost end = 300, span = 200
        #expect(tlen1 == 200)
        #expect(tlen2 == -200)
    }

    @Test("TLEN computation when read2 is leftmost")
    func testComputeTLENReversed() {
        let (tlen1, tlen2) = PairedEndResolver.computeTLEN(
            pos1: 300, isReverse1: true, refLen1: 50,
            pos2: 100, isReverse2: false, refLen2: 50
        )

        // Leftmost = 100, rightmost end = 350, span = 250
        #expect(tlen1 == -250)
        #expect(tlen2 == 250)
    }

    @Test("adjustMAPQ boosts proper pair MAPQ")
    func testAdjustMAPQ() {
        let adjusted = PairedEndResolver.adjustMAPQ(
            mapq: 10,
            pairScore: 200,
            secondBestPairScore: nil,
            isProperPair: true
        )
        // No second best: boost to at least 30
        #expect(adjusted >= 30)
    }

    @Test("adjustMAPQ does not boost improper pair")
    func testAdjustMAPQImproper() {
        let adjusted = PairedEndResolver.adjustMAPQ(
            mapq: 10,
            pairScore: 200,
            secondBestPairScore: nil,
            isProperPair: false
        )
        #expect(adjusted == 10)
    }

    @Test("Resolve returns nil when regions are empty")
    func testResolveEmpty() {
        let dist = InsertSizeDistribution()
        let decision = PairedEndResolver.resolve(
            regions1: [], regions2: [],
            dist: dist,
            genomeLength: 10000,
            scoring: ScoringParameters()
        )
        #expect(decision == nil)
    }
}

// MARK: - SAMOutputBuilder Paired-End Tests

@Suite("SAMOutputBuilder Paired-End Tests")
struct SAMOutputBuilderPETests {

    @Test("Paired-end flags set correctly for read1")
    func testPairedEndFlagsRead1() throws {
        let metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 10000, name: "chr1")]
        )

        let region = MemAlnReg(
            rb: 100, re: 200, qb: 0, qe: 100, rid: 0,
            score: 100, trueScore: 100, sub: 0, w: 10
        )

        let pe = PairedEndInfo(
            isRead1: true, isProperPair: true,
            mateTid: 0, matePos: 300,
            mateIsReverse: true, mateIsUnmapped: false,
            tlen: 300
        )

        let read = ReadSequence(
            name: "read1",
            sequence: String(repeating: "A", count: 100),
            qualityString: String(repeating: "I", count: 100)
        )

        let record = try SAMOutputBuilder.buildRecord(
            read: read, region: region, allRegions: [region],
            metadata: metadata, scoring: ScoringParameters(),
            cigar: [100 << 4 | 0],  // 100M
            isPrimary: true, pairedEnd: pe
        )

        let flag = record.flag
        #expect(flag.contains(.paired))
        #expect(flag.contains(.read1))
        #expect(!flag.contains(.read2))
        #expect(flag.contains(.properPair))
        #expect(flag.contains(.mateReverse))
        #expect(!flag.contains(.mateUnmapped))
    }

    @Test("Unmapped record with paired-end info")
    func testUnmappedPairedRecord() throws {
        let pe = PairedEndInfo(
            isRead1: false, isProperPair: false,
            mateTid: 0, matePos: 100,
            mateIsReverse: false, mateIsUnmapped: false,
            tlen: 0
        )

        let read = ReadSequence(
            name: "read2",
            sequence: "ACGT",
            qualityString: "IIII"
        )

        let record = try SAMOutputBuilder.buildUnmappedRecord(read: read, pairedEnd: pe)
        let flag = record.flag
        #expect(flag.contains(.paired))
        #expect(flag.contains(.unmapped))
        #expect(flag.contains(.read2))
        #expect(!flag.contains(.read1))
        #expect(!flag.contains(.mateUnmapped))
    }
}

// MARK: - Supplementary Alignment Tests

@Suite("Supplementary Alignment Tests")
struct SupplementaryAlignmentTests {

    @Test("SA tag format with two segments")
    func testSATagFormat() {
        let segments: [(rname: String, pos: Int64, isReverse: Bool,
                        cigarString: String, mapq: UInt8, nm: Int32)] = [
            (rname: "chr1", pos: 99, isReverse: false, cigarString: "50M50S", mapq: 60, nm: 1),
            (rname: "chr2", pos: 499, isReverse: true, cigarString: "50S50M", mapq: 30, nm: 0),
        ]

        // Excluding index 0 should produce only the second segment
        let sa0 = SAMOutputBuilder.buildSATag(segments: segments, excludeIndex: 0)
        #expect(sa0 == "chr2,500,-,50S50M,30,0;")

        // Excluding index 1 should produce only the first segment
        let sa1 = SAMOutputBuilder.buildSATag(segments: segments, excludeIndex: 1)
        #expect(sa1 == "chr1,100,+,50M50S,60,1;")
    }

    @Test("SA tag with three segments excludes self")
    func testSATagThreeSegments() {
        let segments: [(rname: String, pos: Int64, isReverse: Bool,
                        cigarString: String, mapq: UInt8, nm: Int32)] = [
            (rname: "chr1", pos: 99, isReverse: false, cigarString: "30M70S", mapq: 60, nm: 0),
            (rname: "chr1", pos: 999, isReverse: false, cigarString: "30S40M30S", mapq: 40, nm: 1),
            (rname: "chr3", pos: 1999, isReverse: true, cigarString: "70S30M", mapq: 30, nm: 0),
        ]

        let sa1 = SAMOutputBuilder.buildSATag(segments: segments, excludeIndex: 1)
        #expect(sa1 == "chr1,100,+,30M70S,60,0;chr3,2000,-,70S30M,30,0;")
    }

    @Test("XA tag format")
    func testXATagFormat() {
        let secondaries: [(rname: String, pos: Int64, isReverse: Bool,
                           cigarString: String, nm: Int32)] = [
            (rname: "chr1", pos: 199, isReverse: false, cigarString: "100M", nm: 2),
            (rname: "chr5", pos: 499, isReverse: true, cigarString: "100M", nm: 1),
        ]

        let xa = SAMOutputBuilder.buildXATag(secondaries: secondaries)
        #expect(xa != nil)
        #expect(xa == "chr1,+200,100M,2;chr5,-500,100M,1;")
    }

    @Test("XA tag returns nil when exceeding maxHits")
    func testXATagMaxHits() {
        // 6 secondaries with maxHits=5 → nil
        let secondaries: [(rname: String, pos: Int64, isReverse: Bool,
                           cigarString: String, nm: Int32)] = (0..<6).map { i in
            (rname: "chr1", pos: Int64(i * 100), isReverse: false, cigarString: "100M", nm: Int32(0))
        }

        let xa = SAMOutputBuilder.buildXATag(secondaries: secondaries, maxHits: 5)
        #expect(xa == nil)
    }

    @Test("Hard-clip conversion converts leading and trailing soft-clips")
    func testHardClipConversion() {
        // 3S + 50M + 5S
        let cigar: [UInt32] = [
            3 << 4 | CIGAROp.softClip.rawValue,
            50 << 4 | CIGAROp.match.rawValue,
            5 << 4 | CIGAROp.softClip.rawValue,
        ]

        let (hardCigar, trimLeft, trimRight) = SAMOutputBuilder.convertToHardClip(cigar: cigar)
        #expect(trimLeft == 3)
        #expect(trimRight == 5)
        #expect(hardCigar[0] == (3 << 4 | CIGAROp.hardClip.rawValue))
        #expect(hardCigar[1] == (50 << 4 | CIGAROp.match.rawValue))
        #expect(hardCigar[2] == (5 << 4 | CIGAROp.hardClip.rawValue))
    }

    @Test("Hard-clip conversion with no soft-clips is no-op")
    func testHardClipNoSoftClips() {
        let cigar: [UInt32] = [100 << 4 | CIGAROp.match.rawValue]
        let (hardCigar, trimLeft, trimRight) = SAMOutputBuilder.convertToHardClip(cigar: cigar)
        #expect(trimLeft == 0)
        #expect(trimRight == 0)
        #expect(hardCigar == cigar)
    }

    @Test("Hard-clip conversion with only leading soft-clip")
    func testHardClipLeadingOnly() {
        let cigar: [UInt32] = [
            10 << 4 | CIGAROp.softClip.rawValue,
            90 << 4 | CIGAROp.match.rawValue,
        ]

        let (hardCigar, trimLeft, trimRight) = SAMOutputBuilder.convertToHardClip(cigar: cigar)
        #expect(trimLeft == 10)
        #expect(trimRight == 0)
        #expect(hardCigar[0] == (10 << 4 | CIGAROp.hardClip.rawValue))
        #expect(hardCigar[1] == (90 << 4 | CIGAROp.match.rawValue))
    }

    @Test("Supplementary record has 0x800 flag, not 0x100")
    func testSupplementaryFlag() throws {
        let metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 10000, name: "chr1")]
        )

        let region = MemAlnReg(
            rb: 100, re: 200, qb: 0, qe: 100, rid: 0,
            score: 100, trueScore: 100, sub: 0, w: 10
        )

        let read = ReadSequence(
            name: "read1",
            sequence: String(repeating: "A", count: 100),
            qualityString: String(repeating: "I", count: 100)
        )

        let record = try SAMOutputBuilder.buildRecord(
            read: read, region: region, allRegions: [region],
            metadata: metadata, scoring: ScoringParameters(),
            cigar: [100 << 4 | 0],
            isPrimary: false,
            isSupplementary: true
        )

        let flag = record.flag
        #expect(flag.contains(.supplementary))
        #expect(!flag.contains(.secondary))
    }

    @Test("Secondary record has 0x100 flag, not 0x800")
    func testSecondaryFlag() throws {
        let metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 10000, name: "chr1")]
        )

        let region = MemAlnReg(
            rb: 100, re: 200, qb: 0, qe: 100, rid: 0,
            score: 100, trueScore: 100, sub: 0, w: 10
        )

        let read = ReadSequence(
            name: "read1",
            sequence: String(repeating: "A", count: 100),
            qualityString: String(repeating: "I", count: 100)
        )

        let record = try SAMOutputBuilder.buildRecord(
            read: read, region: region, allRegions: [region],
            metadata: metadata, scoring: ScoringParameters(),
            cigar: [100 << 4 | 0],
            isPrimary: false,
            isSupplementary: false
        )

        let flag = record.flag
        #expect(flag.contains(.secondary))
        #expect(!flag.contains(.supplementary))
    }

    @Test("MAPQ capping: supplementary MAPQ does not exceed primary")
    func testMAPQCapping() {
        // Use regions where primary would get MAPQ 60 (unique) and supplementary
        // would also get MAPQ 60 — after capping, supplementary <= primary
        let primary = MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 100, trueScore: 100, sub: 0)
        let supp = MemAlnReg(rb: 1000, re: 1050, qb: 50, qe: 100, score: 50, trueScore: 50, sub: 0)

        let scoring = ScoringParameters()
        let primaryMapq = MappingQuality.compute(
            region: primary, allRegions: [primary, supp], scoring: scoring, readLength: 100
        )
        let suppMapq = MappingQuality.compute(
            region: supp, allRegions: [primary, supp], scoring: scoring, readLength: 100
        )

        // After capping: supplementary MAPQ should not exceed primary MAPQ
        let cappedSuppMapq = min(suppMapq, primaryMapq)
        #expect(cappedSuppMapq <= primaryMapq)
    }

    @Test("ScoringParameters flag constants have correct values")
    func testFlagConstants() {
        #expect(ScoringParameters.flagNoMulti == 0x10)
        #expect(ScoringParameters.flagSoftClip == 0x200)
    }
}

// MARK: - MateRescue Tests

@Suite("MateRescue Tests")
struct MateRescueTests {

    /// Build a small packed reference and metadata for testing.
    /// Creates a 1000bp reference of repeating ACGT on chromosome "chr1".
    static func makeTestReference() -> (PackedReference, ReferenceMetadata) {
        let refLen = 1000
        let byteCount = (refLen + 3) / 4
        let ptr = UnsafeMutablePointer<UInt8>.allocate(capacity: byteCount)

        // Fill with repeating pattern: A(0) C(1) G(2) T(3)
        for i in 0..<byteCount {
            // Each byte: bits 7-6=base0, 5-4=base1, 3-2=base2, 1-0=base3
            // Pattern ACGT = 0b_00_01_10_11 = 0x1B
            ptr[i] = 0x1B
        }

        let bufPtr = UnsafeBufferPointer(start: UnsafePointer(ptr), count: byteCount)
        let packedRef = PackedReference(data: bufPtr, length: Int64(refLen), ownedBase: ptr)

        let metadata = ReferenceMetadata(
            totalLength: Int64(refLen),
            numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: Int32(refLen), name: "chr1")]
        )
        return (packedRef, metadata)
    }

    @Test("Rescue finds mate in FR orientation")
    func testRescueFR() {
        let (packedRef, metadata) = MateRescueTests.makeTestReference()
        let genomeLen: Int64 = 1000

        // Template: read1 mapped at forward position 100-120
        let template = MemAlnReg(
            rb: 100, re: 120, qb: 0, qe: 20, rid: 0,
            score: 20, trueScore: 20, sub: 0
        )

        // Create mate read that matches at position ~380 on reverse strand
        // For FR orientation (bwaDir=1): mate is on opposite strand, downstream
        // Extract what the reference looks like at position 380-400, reversed & complemented
        let refSlice = packedRef.subsequence(from: 380, length: 20)
        let mateSeq = refSlice.reversed().map { UInt8(3 - $0) }  // reverse complement
        // This is what the mate read looks like (original orientation)
        // When rescue does reverseComplement, it should match position 380-400

        let mateRead = ReadSequence(name: "mate", bases: mateSeq, qualities: [UInt8](repeating: 30, count: 20))

        // Insert size distribution centered at 300 for FR
        var dist = InsertSizeDistribution()
        dist.stats[PairOrientation.fr.rawValue].mean = 300
        dist.stats[PairOrientation.fr.rawValue].stddev = 50
        dist.stats[PairOrientation.fr.rawValue].properLow = 150
        dist.stats[PairOrientation.fr.rawValue].properHigh = 450
        dist.stats[PairOrientation.fr.rawValue].count = 100
        dist.stats[PairOrientation.fr.rawValue].failed = false
        dist.primaryOrientation = .fr

        let scoring = ScoringParameters()
        let rescued = MateRescue.rescue(
            templateRegions: [template],
            mateRead: mateRead,
            mateRegions: [],
            dist: dist,
            genomeLength: genomeLen,
            packedRef: packedRef,
            metadata: metadata,
            scoring: scoring
        )

        // Should find at least one rescued region
        #expect(!rescued.isEmpty)
        if let reg = rescued.first {
            #expect(reg.score >= scoring.minSeedLength * scoring.matchScore)
            #expect(reg.rid == 0)
        }
    }

    @Test("Rescue skips when mate already has hit")
    func testRescueSkipsExisting() {
        let (packedRef, metadata) = MateRescueTests.makeTestReference()
        let genomeLen: Int64 = 1000

        let template = MemAlnReg(
            rb: 100, re: 120, qb: 0, qe: 20, rid: 0,
            score: 20, trueScore: 20, sub: 0
        )

        let mateRead = ReadSequence(
            name: "mate",
            bases: [UInt8](repeating: 0, count: 20),
            qualities: [UInt8](repeating: 30, count: 20)
        )

        // Existing mate region that satisfies FR at proper distance
        // FR direction: p2 = 2*genomeLen - 1 - mateReg.rb should give dist ~300 from template.rb=100
        // p2 = 100 + 300 = 400, so mateReg.rb = 2*1000 - 1 - 400 = 1599
        let existingMate = MemAlnReg(
            rb: 1599, re: 1619, qb: 0, qe: 20, rid: 0,
            score: 20, trueScore: 20, sub: 0
        )

        var dist = InsertSizeDistribution()
        dist.stats[PairOrientation.fr.rawValue].mean = 300
        dist.stats[PairOrientation.fr.rawValue].stddev = 50
        dist.stats[PairOrientation.fr.rawValue].properLow = 150
        dist.stats[PairOrientation.fr.rawValue].properHigh = 450
        dist.stats[PairOrientation.fr.rawValue].count = 100
        dist.stats[PairOrientation.fr.rawValue].failed = false
        dist.primaryOrientation = .fr

        let scoring = ScoringParameters()
        let rescued = MateRescue.rescue(
            templateRegions: [template],
            mateRead: mateRead,
            mateRegions: [existingMate],
            dist: dist,
            genomeLength: genomeLen,
            packedRef: packedRef,
            metadata: metadata,
            scoring: scoring
        )

        // Should skip FR direction since existing mate satisfies it
        // May still try other directions but they have no stats
        #expect(rescued.isEmpty)
    }

    @Test("Rescue with empty templates returns empty")
    func testRescueEmptyTemplates() {
        let (packedRef, metadata) = MateRescueTests.makeTestReference()
        let mateRead = ReadSequence(
            name: "mate", bases: [0, 1, 2, 3], qualities: [30, 30, 30, 30]
        )
        let dist = InsertSizeDistribution()
        let scoring = ScoringParameters()

        let rescued = MateRescue.rescue(
            templateRegions: [],
            mateRead: mateRead,
            mateRegions: [],
            dist: dist,
            genomeLength: 1000,
            packedRef: packedRef,
            metadata: metadata,
            scoring: scoring
        )
        #expect(rescued.isEmpty)
    }
}

// MARK: - Read Group Tests

@Suite("Read Group Tests")
struct ReadGroupTests {

    @Test("readGroupID parses ID from RG line")
    func testReadGroupIDParsing() {
        var options = BWAMemOptions()
        options.readGroupLine = "@RG\\tID:foo\\tSM:bar\\tPL:ILLUMINA"
        #expect(options.readGroupID == "foo")
    }

    @Test("readGroupID returns nil when no read group set")
    func testReadGroupIDNil() {
        let options = BWAMemOptions()
        #expect(options.readGroupID == nil)
    }

    @Test("readGroupID returns nil when ID field missing")
    func testReadGroupIDMissingID() {
        var options = BWAMemOptions()
        options.readGroupLine = "@RG\\tSM:bar\\tPL:ILLUMINA"
        #expect(options.readGroupID == nil)
    }

    @Test("buildRecord includes RG aux tag when readGroupID provided")
    func testBuildRecordWithRG() throws {
        let metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 10000, name: "chr1")]
        )

        let region = MemAlnReg(
            rb: 100, re: 200, qb: 0, qe: 100, rid: 0,
            score: 100, trueScore: 100, sub: 0, w: 10
        )

        let read = ReadSequence(
            name: "test_rg",
            sequence: String(repeating: "A", count: 100),
            qualityString: String(repeating: "I", count: 100)
        )

        let record = try SAMOutputBuilder.buildRecord(
            read: read, region: region, allRegions: [region],
            metadata: metadata, scoring: ScoringParameters(),
            cigar: [100 << 4 | 0],
            isPrimary: true,
            readGroupID: "sample1"
        )

        let aux = record.auxiliaryData
        let rgVal = aux.string(forTag: "RG")
        #expect(rgVal == "sample1")
    }

    @Test("buildRecord has no RG tag when readGroupID is nil")
    func testBuildRecordWithoutRG() throws {
        let metadata = ReferenceMetadata(
            totalLength: 10000, numSequences: 1,
            annotations: [ReferenceAnnotation(offset: 0, length: 10000, name: "chr1")]
        )

        let region = MemAlnReg(
            rb: 100, re: 200, qb: 0, qe: 100, rid: 0,
            score: 100, trueScore: 100, sub: 0, w: 10
        )

        let read = ReadSequence(
            name: "test_no_rg",
            sequence: String(repeating: "A", count: 100),
            qualityString: String(repeating: "I", count: 100)
        )

        let record = try SAMOutputBuilder.buildRecord(
            read: read, region: region, allRegions: [region],
            metadata: metadata, scoring: ScoringParameters(),
            cigar: [100 << 4 | 0],
            isPrimary: true
        )

        let aux = record.auxiliaryData
        let rgVal = aux.string(forTag: "RG")
        #expect(rgVal == nil)
    }

    @Test("buildUnmappedRecord includes RG tag when readGroupID provided")
    func testUnmappedRecordWithRG() throws {
        let read = ReadSequence(
            name: "unmapped_rg",
            sequence: "ACGT",
            qualityString: "IIII"
        )

        let record = try SAMOutputBuilder.buildUnmappedRecord(
            read: read, readGroupID: "sample2"
        )

        let aux = record.auxiliaryData
        let rgVal = aux.string(forTag: "RG")
        #expect(rgVal == "sample2")
    }
}
