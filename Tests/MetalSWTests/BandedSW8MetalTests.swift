#if canImport(Metal)
import Testing
@testable import MetalSW
@testable import Alignment
@testable import BWACore

@Suite("BandedSW8 Metal Tests")
struct BandedSW8MetalTests {

    @Test("Metal engine initializes")
    func testEngineInit() {
        let engine = MetalSWEngine.shared
        #expect(engine != nil, "MetalSWEngine should initialize on macOS with Metal")
    }

    @Test("Perfect match — GPU matches CPU")
    func testPerfectMatch() throws {
        guard let engine = MetalSWEngine.shared else {
            return  // Skip on systems without Metal
        }

        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]  // ACGTACGTAC
        let target: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1]

        // CPU result
        let cpuResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: 0)
            }
        }

        // GPU result
        let task = BandedSWTask(query: query, target: target, h0: 0, w: 100, scoring: scoring)
        let gpuResults = BandedSWDispatcher.dispatchBatch8(tasks: [task], engine: engine)

        #expect(gpuResults.count == 1)
        if let cpuR = cpuResult, let gpuR = gpuResults[0] {
            #expect(gpuR.score == cpuR.score, "Score mismatch: GPU=\(gpuR.score) CPU=\(cpuR.score)")
            #expect(gpuR.queryEnd == cpuR.queryEnd)
            #expect(gpuR.targetEnd == cpuR.targetEnd)
            #expect(gpuR.globalScore == cpuR.globalScore)
            #expect(gpuR.globalTargetEnd == cpuR.globalTargetEnd)
        } else {
            #expect(cpuResult != nil, "CPU result should not be nil for perfect match")
            #expect(gpuResults[0] != nil, "GPU result should not be nil for perfect match")
        }
    }

    @Test("Mismatch — GPU matches CPU")
    func testMismatch() throws {
        guard let engine = MetalSWEngine.shared else { return }

        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3]   // ACGTACGT
        let target: [UInt8] = [0, 1, 0, 3, 0, 1, 0, 3]   // ACATACGT -> 2 mismatches

        let cpuResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: 0)
            }
        }

        let task = BandedSWTask(query: query, target: target, h0: 0, w: 100, scoring: scoring)
        let gpuResults = BandedSWDispatcher.dispatchBatch8(tasks: [task], engine: engine)

        #expect(gpuResults.count == 1)
        if let cpuR = cpuResult, let gpuR = gpuResults[0] {
            #expect(gpuR.score == cpuR.score)
            #expect(gpuR.queryEnd == cpuR.queryEnd)
            #expect(gpuR.targetEnd == cpuR.targetEnd)
        }
    }

    @Test("With h0 > 0 — GPU matches CPU")
    func testWithH0() throws {
        guard let engine = MetalSWEngine.shared else { return }

        let scoring = ScoringParameters()
        let query: [UInt8] = [0, 1, 2, 3, 0, 1]
        let target: [UInt8] = [0, 1, 2, 3, 0, 1]

        let cpuResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: 10)
            }
        }

        let task = BandedSWTask(query: query, target: target, h0: 10, w: 100, scoring: scoring)
        let gpuResults = BandedSWDispatcher.dispatchBatch8(tasks: [task], engine: engine)

        #expect(gpuResults.count == 1)
        if let cpuR = cpuResult, let gpuR = gpuResults[0] {
            #expect(gpuR.score == cpuR.score)
            #expect(gpuR.queryEnd == cpuR.queryEnd)
            #expect(gpuR.targetEnd == cpuR.targetEnd)
        }
    }

    @Test("Batch of multiple tasks — all match CPU")
    func testBatch() throws {
        guard let engine = MetalSWEngine.shared else { return }

        let scoring = ScoringParameters()
        let testCases: [(query: [UInt8], target: [UInt8], h0: Int32)] = [
            ([0, 1, 2, 3], [0, 1, 2, 3], 0),
            ([0, 1, 2, 3, 0, 1, 2, 3], [0, 1, 2, 3, 0, 1, 2, 3], 5),
            ([0, 0, 0, 0, 0], [0, 0, 0, 0, 0], 0),
            ([0, 1, 2, 3, 0, 1], [3, 2, 1, 0, 3, 2], 0),  // all mismatches
        ]

        var tasks: [BandedSWTask] = []
        var cpuResults: [SWResult?] = []

        for tc in testCases {
            tasks.append(BandedSWTask(query: tc.query, target: tc.target, h0: tc.h0, w: 100, scoring: scoring))
            let cpu = tc.query.withUnsafeBufferPointer { qBuf in
                tc.target.withUnsafeBufferPointer { tBuf in
                    BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: tc.h0)
                }
            }
            cpuResults.append(cpu)
        }

        let gpuResults = BandedSWDispatcher.dispatchBatch8(tasks: tasks, engine: engine)

        #expect(gpuResults.count == testCases.count)
        for i in 0..<testCases.count {
            if let cpuR = cpuResults[i], let gpuR = gpuResults[i] {
                #expect(gpuR.score == cpuR.score, "Task \(i): score GPU=\(gpuR.score) CPU=\(cpuR.score)")
                #expect(gpuR.queryEnd == cpuR.queryEnd, "Task \(i): queryEnd mismatch")
                #expect(gpuR.targetEnd == cpuR.targetEnd, "Task \(i): targetEnd mismatch")
            } else {
                // Both should be nil (overflow) or both non-nil
                #expect(cpuResults[i] == nil && gpuResults[i] == nil,
                        "Task \(i): nil mismatch cpu=\(cpuResults[i] != nil) gpu=\(gpuResults[i] != nil)")
            }
        }
    }

    @Test("150bp read — GPU matches CPU")
    func test150bpRead() throws {
        guard let engine = MetalSWEngine.shared else { return }

        let scoring = ScoringParameters()
        // Simulate a realistic 150bp read
        var query: [UInt8] = []
        var target: [UInt8] = []
        for i in 0..<150 {
            query.append(UInt8(i % 4))
            target.append(UInt8(i % 4))
        }
        // Introduce a few mismatches
        target[50] = (target[50] + 1) % 4
        target[100] = (target[100] + 2) % 4

        let cpuResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: 0)
            }
        }

        let task = BandedSWTask(query: query, target: target, h0: 0, w: 100, scoring: scoring)
        let gpuResults = BandedSWDispatcher.dispatchBatch8(tasks: [task], engine: engine)

        if let cpuR = cpuResult, let gpuR = gpuResults[0] {
            #expect(gpuR.score == cpuR.score)
            #expect(gpuR.queryEnd == cpuR.queryEnd)
            #expect(gpuR.targetEnd == cpuR.targetEnd)
            #expect(gpuR.globalScore == cpuR.globalScore)
            // maxOff may differ: scalar GPU iterates linearly while CPU uses Farrar striped
            // order, so when multiple positions tie for max score, different j is found first
        }
    }

    @Test("Overflow signals nil (fall back to 16-bit)")
    func testOverflow() throws {
        guard let engine = MetalSWEngine.shared else { return }

        var scoring = ScoringParameters()
        scoring.matchScore = 10  // high match score to force overflow
        let query: [UInt8] = Array(repeating: 0, count: 50)
        let target: [UInt8] = Array(repeating: 0, count: 50)

        let cpuResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW8.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: 0)
            }
        }

        let task = BandedSWTask(query: query, target: target, h0: 0, w: 100, scoring: scoring)
        let gpuResults = BandedSWDispatcher.dispatchBatch8(tasks: [task], engine: engine)

        // Both should overflow (return nil)
        #expect(cpuResult == nil, "CPU should detect overflow")
        #expect(gpuResults[0] == nil, "GPU should detect overflow")
    }

    @Test("8-bit overflow falls back to 16-bit GPU")
    func testOverflowFallbackTo16bit() throws {
        guard let engine = MetalSWEngine.shared else { return }
        guard engine.bandedSW16Pipeline != nil else { return }

        var scoring = ScoringParameters()
        scoring.matchScore = 10
        let query: [UInt8] = Array(repeating: 0, count: 50)
        let target: [UInt8] = Array(repeating: 0, count: 50)

        // CPU 16-bit result
        let cpuResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: 0)
            }
        }

        // GPU 16-bit result
        let task = BandedSWTask(query: query, target: target, h0: 0, w: 100, scoring: scoring)
        let gpuResults = BandedSWDispatcher.dispatchBatch16(tasks: [task], engine: engine)

        #expect(gpuResults.count == 1)
        let gpuR = gpuResults[0]
        #expect(gpuR.score == cpuResult.score, "Score: GPU=\(gpuR.score) CPU=\(cpuResult.score)")
        #expect(gpuR.queryEnd == cpuResult.queryEnd)
        #expect(gpuR.targetEnd == cpuResult.targetEnd)
    }

    @Test("16-bit with mismatches — GPU matches CPU")
    func test16bitMismatch() throws {
        guard let engine = MetalSWEngine.shared else { return }
        guard engine.bandedSW16Pipeline != nil else { return }

        var scoring = ScoringParameters()
        scoring.matchScore = 8
        let query: [UInt8] = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3,
                              0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
        var target = query
        target[5] = 3  // mismatch
        target[15] = 0  // mismatch

        let cpuResult = query.withUnsafeBufferPointer { qBuf in
            target.withUnsafeBufferPointer { tBuf in
                BandedSW16.align(query: qBuf, target: tBuf, scoring: scoring, w: 100, h0: 0)
            }
        }

        let task = BandedSWTask(query: query, target: target, h0: 0, w: 100, scoring: scoring)
        let gpuResults = BandedSWDispatcher.dispatchBatch16(tasks: [task], engine: engine)

        #expect(gpuResults.count == 1)
        let gpuR = gpuResults[0]
        #expect(gpuR.score == cpuResult.score)
        #expect(gpuR.queryEnd == cpuResult.queryEnd)
        #expect(gpuR.targetEnd == cpuResult.targetEnd)
    }
}
#endif
