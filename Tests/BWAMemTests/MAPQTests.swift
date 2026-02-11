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

        // Region with sub=0 (current behavior for unique mappers)
        let regionNoSub = MemAlnReg(
            rb: 0, re: 100, qb: 0, qe: 100, score: 100, sub: 0
        )

        // Same region but with sub=80 (what it should be if a secondary existed)
        let regionWithSub = MemAlnReg(
            rb: 0, re: 100, qb: 0, qe: 100, score: 100, sub: 80
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

        // With sub=0: scoreDiff=100, ratio=1.0, identity=0.667 → capped at 60
        // With sub=80: scoreDiff=20, ratio=0.2, identity=0.667 → ~20
        #expect(mapqWithSub < mapqNoSub,
                "sub=0 inflates MAPQ: \(mapqNoSub) vs \(mapqWithSub) with sub=80")
        // sub=0 should give very high MAPQ (capped at 60)
        #expect(mapqNoSub >= 50, "Unique mapper with sub=0 gets near-max MAPQ")
        // sub=80 should give noticeably lower MAPQ
        #expect(mapqWithSub < 40, "Sub=80 should substantially reduce MAPQ")
    }
}
