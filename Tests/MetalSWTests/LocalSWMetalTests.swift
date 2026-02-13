#if canImport(Metal)
import Testing
@testable import MetalSW
@testable import Alignment
@testable import BWACore

@Suite("LocalSW Metal Tests")
struct LocalSWMetalTests {

    @Test("Perfect match — GPU matches CPU")
    func testPerfectMatch() async throws {
        guard let engine = MetalSWEngine.shared else { return }
        guard engine.localSWPipeline != nil else { return }

        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]

        let cpuResult = LocalSWAligner.align(query: query, target: target, scoring: scoring)

        let task = LocalSWTask(query: query, target: target, scoring: scoring)
        let gpuResults = await LocalSWDispatcher.dispatchBatch(tasks: [task], engine: engine)

        #expect(gpuResults.count == 1)
        if let cpuR = cpuResult, let gpuR = gpuResults[0] {
            #expect(gpuR.score == cpuR.score, "Score: GPU=\(gpuR.score) CPU=\(cpuR.score)")
            #expect(gpuR.queryBegin == cpuR.queryBegin)
            #expect(gpuR.queryEnd == cpuR.queryEnd)
            #expect(gpuR.targetBegin == cpuR.targetBegin)
            #expect(gpuR.targetEnd == cpuR.targetEnd)
        }
    }

    @Test("Subsequence match — GPU matches CPU")
    func testSubsequenceMatch() async throws {
        guard let engine = MetalSWEngine.shared else { return }
        guard engine.localSWPipeline != nil else { return }

        let scoring = ScoringParameters()
        // Query is a substring of target with flanking mismatches
        let query: [UInt8] =  [0, 1, 2, 3, 0, 1]
        let target: [UInt8] = [3, 3, 0, 1, 2, 3, 0, 1, 3, 3]

        let cpuResult = LocalSWAligner.align(query: query, target: target, scoring: scoring)

        let task = LocalSWTask(query: query, target: target, scoring: scoring)
        let gpuResults = await LocalSWDispatcher.dispatchBatch(tasks: [task], engine: engine)

        #expect(gpuResults.count == 1)
        if let cpuR = cpuResult, let gpuR = gpuResults[0] {
            #expect(gpuR.score == cpuR.score)
            #expect(gpuR.queryBegin == cpuR.queryBegin)
            #expect(gpuR.queryEnd == cpuR.queryEnd)
            #expect(gpuR.targetBegin == cpuR.targetBegin)
            #expect(gpuR.targetEnd == cpuR.targetEnd)
        }
    }

    @Test("Batch of multiple tasks — all match CPU")
    func testBatch() async throws {
        guard let engine = MetalSWEngine.shared else { return }
        guard engine.localSWPipeline != nil else { return }

        let scoring = ScoringParameters()
        let testCases: [(query: [UInt8], target: [UInt8])] = [
            ([0, 1, 2, 3], [0, 1, 2, 3]),
            ([0, 1, 2, 3, 0, 1, 2, 3], [0, 1, 2, 3, 0, 1, 2, 3]),
            ([0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
            ([0, 1, 2, 3], [3, 2, 1, 0]),  // all mismatches
        ]

        var tasks: [LocalSWTask] = []
        var cpuResults: [Alignment.LocalSWResult?] = []

        for tc in testCases {
            tasks.append(LocalSWTask(query: tc.query, target: tc.target, scoring: scoring))
            cpuResults.append(LocalSWAligner.align(query: tc.query, target: tc.target, scoring: scoring))
        }

        let gpuResults = await LocalSWDispatcher.dispatchBatch(tasks: tasks, engine: engine)

        for i in 0..<testCases.count {
            if let cpuR = cpuResults[i], let gpuR = gpuResults[i] {
                #expect(gpuR.score == cpuR.score, "Task \(i): score GPU=\(gpuR.score) CPU=\(cpuR.score)")
                // Coordinates may differ on tie-breaking (multiple cells with same max score).
                // Only require coordinate match for non-degenerate alignments (score > 1).
                if cpuR.score > 1 {
                    #expect(gpuR.queryBegin == cpuR.queryBegin, "Task \(i): queryBegin mismatch")
                    #expect(gpuR.queryEnd == cpuR.queryEnd, "Task \(i): queryEnd mismatch")
                }
            } else if cpuResults[i] == nil && gpuResults[i] == nil {
                // Both nil (no alignment or overflow) — OK
            } else if cpuResults[i] != nil && gpuResults[i] == nil {
                // GPU returned nil but CPU found an alignment — check if it's a zero score
                #expect(cpuResults[i]!.score == 0, "Task \(i): GPU nil but CPU had score \(cpuResults[i]!.score)")
            }
        }
    }

    @Test("150bp mate rescue scenario — GPU matches CPU")
    func test150bpMateRescue() async throws {
        guard let engine = MetalSWEngine.shared else { return }
        guard engine.localSWPipeline != nil else { return }

        let scoring = ScoringParameters()
        var query: [UInt8] = []
        var target: [UInt8] = []
        // Simulate mate rescue: query is 150bp read, target is ~500bp region
        for i in 0..<150 {
            query.append(UInt8(i % 4))
        }
        for i in 0..<500 {
            target.append(UInt8(i % 4))
        }
        // Mismatches
        target[200] = (target[200] + 1) % 4
        target[250] = (target[250] + 2) % 4

        let cpuResult = LocalSWAligner.align(query: query, target: target, scoring: scoring)

        let task = LocalSWTask(query: query, target: target, scoring: scoring)
        let gpuResults = await LocalSWDispatcher.dispatchBatch(tasks: [task], engine: engine)

        if let cpuR = cpuResult, let gpuR = gpuResults[0] {
            #expect(gpuR.score == cpuR.score)
            #expect(gpuR.queryBegin == cpuR.queryBegin)
            #expect(gpuR.queryEnd == cpuR.queryEnd)
            #expect(gpuR.targetBegin == cpuR.targetBegin)
            #expect(gpuR.targetEnd == cpuR.targetEnd)
        }
    }
}
#endif
