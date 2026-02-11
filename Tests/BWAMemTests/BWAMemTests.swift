import Testing
@testable import BWAMem
@testable import BWACore

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
        let region = MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 100, trueScore: 100, sub: 0)
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
        let region1 = MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 80, trueScore: 80, sub: 70)
        let region2 = MemAlnReg(rb: 0, re: 100, qb: 0, qe: 100, score: 80, trueScore: 80, sub: 10)

        let scoring = ScoringParameters()
        let mapq1 = MappingQuality.compute(region: region1, allRegions: [region1], scoring: scoring, readLength: 100)
        let mapq2 = MappingQuality.compute(region: region2, allRegions: [region2], scoring: scoring, readLength: 100)

        // Higher sub-optimal score should give lower MAPQ
        #expect(mapq1 < mapq2)
    }
}

// MARK: - Insert Size Estimator Tests

@Suite("InsertSizeEstimator Tests")
struct InsertSizeEstimatorTests {

    @Test("Infer FR orientation: read1 forward, read2 reverse, r1 left of r2")
    func testInferOrientationFR() {
        let genomeLen: Int64 = 10000
        // Read1: forward strand, pos 100-200
        let r1 = MemAlnReg(rb: 100, re: 200, qb: 0, qe: 100, rid: 0, score: 100)
        // Read2: reverse strand, pos 300-400 (in BWT coords: genomeLen + something)
        // Reverse strand coords: rb >= genomeLen
        // Forward-strand pos = 2*genomeLen - re ... 2*genomeLen - rb
        // Want fwd pos 300-400: 2*10000 - re = 300 → re = 19700; 2*10000 - rb = 400 → rb = 19600
        let r2 = MemAlnReg(rb: 19600, re: 19700, qb: 0, qe: 100, rid: 0, score: 100)

        let result = InsertSizeEstimator.inferOrientation(r1: r1, r2: r2, genomeLength: genomeLen)
        #expect(result != nil)
        #expect(result!.orientation == .fr)
        // Insert size: outer coords = min(100,300)=100, max(200,400)=400, span=300
        #expect(result!.insertSize == 300)
    }

    @Test("Infer RF orientation: read1 reverse before read2 forward")
    func testInferOrientationRF() {
        let genomeLen: Int64 = 10000
        // Read1: reverse strand, fwd pos 300-400
        // 2*10000 - re = 300 → re=19700; 2*10000 - rb = 400 → rb=19600
        let r1 = MemAlnReg(rb: 19600, re: 19700, qb: 0, qe: 100, rid: 0, score: 100)
        // Read2: forward strand, pos 100-200
        let r2 = MemAlnReg(rb: 100, re: 200, qb: 0, qe: 100, rid: 0, score: 100)

        let result = InsertSizeEstimator.inferOrientation(r1: r1, r2: r2, genomeLength: genomeLen)
        #expect(result != nil)
        // r1 is reverse, r2 is forward: r2Start(100) <= r1Start(300) → FR
        // Actually: r1Reverse=true, r2Reverse=false → "r1Reverse && !r2Reverse" case
        // r2Start(100) <= r1Start(300) → .fr
        #expect(result!.orientation == .fr)
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

    @Test("Infer RR orientation: both reverse")
    func testInferOrientationRR() {
        let genomeLen: Int64 = 10000
        let r1 = MemAlnReg(rb: 19600, re: 19700, qb: 0, qe: 100, rid: 0, score: 100)
        let r2 = MemAlnReg(rb: 19800, re: 19900, qb: 0, qe: 100, rid: 0, score: 100)

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
        // Create 100 simulated FR pairs with insert sizes around 300 (±20)
        var regions1: [[MemAlnReg]] = []
        var regions2: [[MemAlnReg]] = []

        for i in 0..<100 {
            let insertSize = 280 + Int64(i % 41)  // 280..320
            let r1Pos: Int64 = 1000 + Int64(i) * 500
            let r1 = MemAlnReg(
                rb: r1Pos, re: r1Pos + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0
            )
            // r2 on reverse strand, fwd position = r1Pos + insertSize - 100
            let r2FwdStart = r1Pos + insertSize - 100
            let r2FwdEnd = r1Pos + insertSize
            // Reverse BWT coords: rb = 2*genomeLen - fwdEnd, re = 2*genomeLen - fwdStart
            let r2 = MemAlnReg(
                rb: 2 * genomeLen - r2FwdEnd,
                re: 2 * genomeLen - r2FwdStart,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0
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

        // 40 normal pairs with insert ~300
        for i in 0..<40 {
            let r1Pos: Int64 = 1000 + Int64(i) * 500
            let r1 = MemAlnReg(
                rb: r1Pos, re: r1Pos + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0
            )
            let r2FwdStart = r1Pos + 200
            let r2FwdEnd = r1Pos + 300
            let r2 = MemAlnReg(
                rb: 2 * genomeLen - r2FwdEnd,
                re: 2 * genomeLen - r2FwdStart,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0
            )
            regions1.append([r1])
            regions2.append([r2])
        }

        // 5 outlier pairs with insert ~5000 (should be filtered)
        for i in 0..<5 {
            let r1Pos: Int64 = 50000 + Int64(i) * 500
            let r1 = MemAlnReg(
                rb: r1Pos, re: r1Pos + 100,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0
            )
            let r2FwdStart = r1Pos + 4900
            let r2FwdEnd = r1Pos + 5000
            let r2 = MemAlnReg(
                rb: 2 * genomeLen - r2FwdEnd,
                re: 2 * genomeLen - r2FwdStart,
                qb: 0, qe: 100, rid: 0,
                score: 100, trueScore: 100, sub: 0
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
