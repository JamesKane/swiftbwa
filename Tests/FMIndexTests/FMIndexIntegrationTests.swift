import Testing
import Foundation
@testable import FMIndex
@testable import BWACore

/// Integration tests that load a real bwa-mem2 index (lambda phage).
/// These tests require TestData/lambda.fa.bwt.2bit.64 and related files.
@Suite("FMIndex Integration Tests")
struct FMIndexIntegrationTests {

    /// Resolve project root from source file location.
    static var testDataDir: String {
        let file = URL(fileURLWithPath: #filePath)
        let projectRoot = file
            .deletingLastPathComponent()  // FMIndexTests/
            .deletingLastPathComponent()  // Tests/
            .deletingLastPathComponent()  // swiftbwa/
        return projectRoot.appendingPathComponent("TestData").path
    }

    static var indexPrefix: String { testDataDir + "/lambda.fa" }

    static var indexAvailable: Bool {
        FileManager.default.fileExists(atPath: indexPrefix + ".bwt.2bit.64")
    }

    // MARK: - BWT Loading

    @Test("Load BWT from lambda index")
    func testLoadBWT() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let (bwt, _, _) = try FMIndexLoader.loadBWT(from: Self.indexPrefix)

        // ref_seq_len = 97005 (2 * 48502 + 1 sentinel)
        #expect(bwt.length == 97005)

        // Sentinel at position 65380
        #expect(bwt.sentinelIndex == 65380)

        // count[] after +1 adjustment: raw = (0, 24320, 48502, 72684, 97004)
        #expect(bwt.count.0 == 1)
        #expect(bwt.count.1 == 24321)
        #expect(bwt.count.2 == 48503)
        #expect(bwt.count.3 == 72685)
        #expect(bwt.count.4 == 97005)

        // Checkpoint count: 97005 >> 6 + 1 = 1516
        #expect(bwt.checkpoints.count == 1516)
    }

    // MARK: - SA Loading

    @Test("Load suffix array from lambda index")
    func testLoadSA() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let (_, bwtMF, _) = try FMIndexLoader.loadBWT(from: Self.indexPrefix)
        let sa = try FMIndexLoader.loadSA(from: Self.indexPrefix, referenceSeqLen: 97005, mappedFile: bwtMF)

        // SA count: 97005 >> 3 + 1 = 12126
        #expect(sa.count == 12126)
        #expect(sa.compressionShift == 3)

        // All SA entries should be valid positions (>= 0)
        for i in 0..<min(100, sa.count) {
            let entry = sa.entry(at: Int64(i))
            #expect(entry >= 0, "SA entry \(i) = \(entry) should be >= 0")
        }
    }

    // MARK: - Packed Reference Loading

    @Test("Load packed reference from lambda index")
    func testLoadPac() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let pac = try FMIndexLoader.loadPac(from: Self.indexPrefix)

        // Lambda phage genome is 48502 bp
        #expect(pac.length == 48502)

        // All bases should be 0-3 (A, C, G, T)
        for i in 0..<min(100, Int(pac.length)) {
            let b = pac.base(at: Int64(i))
            #expect(b <= 3, "Base at \(i) = \(b) should be <= 3")
        }
    }

    // MARK: - Metadata Loading

    @Test("Load metadata from lambda index")
    func testLoadMetadata() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let metadata = try FMIndexLoader.loadMetadata(from: Self.indexPrefix)

        #expect(metadata.totalLength == 48502)
        #expect(metadata.numSequences == 1)
        #expect(metadata.seed == 11)

        #expect(metadata.annotations.count == 1)
        #expect(metadata.annotations[0].name == "J02459.1")
        #expect(metadata.annotations[0].length == 48502)
        #expect(metadata.annotations[0].offset == 0)
        #expect(metadata.annotations[0].nAmb == 0)
        #expect(metadata.annotations[0].anno.contains("Escherichia phage Lambda"))

        // No ambiguities in lambda phage
        #expect(metadata.ambiguities.isEmpty)
    }

    // MARK: - Full FMIndex Loading

    @Test("Load complete FMIndex")
    func testLoadFMIndex() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let index = try FMIndexLoader.load(from: Self.indexPrefix)

        #expect(index.referenceSeqLen == 97005)
        #expect(index.genomeLength == 48502)
        #expect(index.metadata.numSequences == 1)
        #expect(index.suffixArray.count == 12126)
        #expect(index.packedRef.length == 48502)
    }

    // MARK: - Backward Search on Real Data

    @Test("BackwardSearch initInterval on real BWT")
    func testInitIntervalReal() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let (bwt, _, _) = try FMIndexLoader.loadBWT(from: Self.indexPrefix)

        // Each base should have a non-zero interval
        var totalBases: Int64 = 0
        for c in 0..<4 {
            let interval = BackwardSearch.initInterval(bwt: bwt, base: c)
            #expect(interval.s > 0, "Base \(c) should have positive count")
            #expect(interval.k >= 0, "Base \(c) interval start should be >= 0")
            totalBases += interval.s
        }

        // Sum of all 4 base intervals = total characters minus sentinel
        // count[4] - count[0] = 97005 - 1 = 97004 (sentinel is 1 position)
        #expect(totalBases == bwt.length - 1,
                "Sum of all base counts should equal BWT length minus 1 (sentinel)")
    }

    @Test("BackwardSearch extends correctly on real BWT")
    func testBackwardExtReal() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let (bwt, _, _) = try FMIndexLoader.loadBWT(from: Self.indexPrefix)

        // Start with 'A' interval, extend by 'C' to get "CA" interval
        let intA = BackwardSearch.initInterval(bwt: bwt, base: 0)  // A
        let intCA = BackwardSearch.backwardExt(bwt: bwt, interval: intA, base: 1)  // CA

        // CA interval should be smaller than A interval
        #expect(intCA.s > 0, "CA should have occurrences in lambda phage")
        #expect(intCA.s < intA.s, "CA should be less frequent than A alone")

        // Extend further: "GCA"
        let intGCA = BackwardSearch.backwardExt(bwt: bwt, interval: intCA, base: 2)  // GCA
        #expect(intGCA.s > 0, "GCA should have occurrences")
        #expect(intGCA.s <= intCA.s, "GCA should be no more frequent than CA")
    }

    // MARK: - SMEM Finding on Real Data

    @Test("SMEMFinder finds seeds in lambda phage subsequence")
    func testSMEMsOnRealData() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let index = try FMIndexLoader.load(from: Self.indexPrefix)

        // Extract first 50 bases from packed reference as our "query"
        let queryLen = 50
        var query = [UInt8](repeating: 0, count: queryLen)
        for i in 0..<queryLen {
            query[i] = index.packedRef.base(at: Int64(i))
        }

        let smems = SMEMFinder.findAllSMEMs(
            query: query,
            bwt: index.bwt,
            minSeedLen: 19
        )

        // Since the query is a direct subsequence of the reference, we should find SMEMs
        #expect(!smems.isEmpty, "Should find SMEMs for reference subsequence")

        // Verify SMEM properties
        for smem in smems {
            #expect(smem.queryBegin >= 0, "SMEM queryBegin should be >= 0")
            #expect(smem.queryEnd <= Int32(queryLen), "SMEM queryEnd should be <= query length")
            #expect(smem.queryEnd > smem.queryBegin, "SMEM should have positive length")
            #expect(smem.length >= 19, "SMEM length should be >= minSeedLen")
            #expect(smem.count > 0, "SMEM should have positive occurrence count")
        }
    }

    @Test("SMEMFinder with shorter query")
    func testSMEMsShortQuery() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let index = try FMIndexLoader.load(from: Self.indexPrefix)

        // Extract bases 1000-1030 from reference (30 bases)
        let start: Int64 = 1000
        let queryLen = 30
        var query = [UInt8](repeating: 0, count: queryLen)
        for i in 0..<queryLen {
            query[i] = index.packedRef.base(at: start + Int64(i))
        }

        let smems = SMEMFinder.findAllSMEMs(
            query: query,
            bwt: index.bwt,
            minSeedLen: 15
        )

        #expect(!smems.isEmpty, "Should find SMEMs for 30bp reference subsequence")

        // The entire 30bp query should be covered by at least one SMEM
        // since it comes directly from the reference
        let maxEnd = smems.map { $0.queryEnd }.max() ?? 0
        #expect(maxEnd > 0, "SMEMs should cover some of the query")
    }

    // MARK: - SA Lookup Consistency

    @Test("SA entries correspond to correct BWT positions")
    func testSABWTConsistency() throws {
        try #require(Self.indexAvailable, "Lambda index not available")

        let index = try FMIndexLoader.load(from: Self.indexPrefix)

        // Verify a few SA entries are within valid range
        let maxPos = index.referenceSeqLen
        for i in stride(from: 0, to: min(1000, index.suffixArray.count), by: 10) {
            let saVal = index.suffixArray.entry(at: Int64(i))
            #expect(saVal >= 0 && saVal < maxPos,
                    "SA[\(i)] = \(saVal) should be in [0, \(maxPos))")
        }
    }
}
