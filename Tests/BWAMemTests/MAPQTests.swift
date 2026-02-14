import Testing
@testable import BWAMem
@testable import BWACore

@Suite("MAPQ Divergence Tests")
struct MAPQTests {

    @Test("MAPQ dramatically higher when sub is zero vs nonzero")
    func testMAPQWithAndWithoutSub() {
        // Documents that sub=0 (the default for unique mappers) produces
        // dramatically inflated MAPQ compared to having a realistic sub score.
        //
        // Use readLength=150 with score=100 to avoid the "perfect unique hit"
        // shortcut (score == maxScore → 60) and exercise the full MAPQ formula.
        let scoring = ScoringParameters()
        let readLength: Int32 = 150

        // Region with sub=0 (unique mapper: defaults to minSeedLen * matchScore = 19)
        let regionNoSub = MemAlnReg(
            rb: 0, re: 100, qb: 0, qe: 100, score: 100, sub: 0, seedCov: 100
        )

        // Same region but with sub=80 (what it should be if a secondary existed)
        let regionWithSub = MemAlnReg(
            rb: 0, re: 100, qb: 0, qe: 100, score: 100, sub: 80, seedCov: 100
        )

        let mapqNoSub = MappingQuality.compute(
            region: regionNoSub,
            allRegions: [regionNoSub],
            scoring: scoring,
            readLength: readLength
        )

        let mapqWithSub = MappingQuality.compute(
            region: regionWithSub,
            allRegions: [regionWithSub],
            scoring: scoring,
            readLength: readLength
        )

        // bwa-mem2's mapQ_coef_len=50 formula: 6.02 * (score-sub)/a * tmp^2
        // With sub=0 (defaults to 19): 6.02 * 81 * 0.424 = 206 → capped at 60
        // With sub=80: 6.02 * 20 * 0.424 = 51
        #expect(mapqWithSub < mapqNoSub,
                "sub=0 inflates MAPQ: \(mapqNoSub) vs \(mapqWithSub) with sub=80")
        // sub=0 should give very high MAPQ (capped at 60)
        #expect(mapqNoSub >= 50, "Unique mapper with sub=0 gets near-max MAPQ")
        // sub=80 should give lower MAPQ than sub=0
        #expect(mapqWithSub < 55, "Sub=80 should reduce MAPQ below sub=0 level")
    }
}
