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

    @Test("InsertSizeEstimator running statistics")
    func testInsertSizeEstimator() {
        var estimator = InsertSizeEstimator()

        // Add insert sizes around 300 with some variance
        for size in [280, 290, 300, 310, 320] {
            estimator.add(Int64(size))
        }

        #expect(estimator.count == 5)
        #expect(estimator.mean == 300.0)
        #expect(estimator.stddev > 0)
    }

    @Test("InsertSizeEstimator proper pair check")
    func testProperPair() {
        var estimator = InsertSizeEstimator()

        // Need 25+ samples for proper pair checking
        for _ in 0..<30 {
            estimator.add(300)
        }

        #expect(estimator.isProperPair(300) == true)
        #expect(estimator.isProperPair(5000) == false)  // Way outside range
    }
}
