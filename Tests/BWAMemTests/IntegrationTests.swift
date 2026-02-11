import Foundation
import Testing
@testable import BWAMem
@testable import BWACore
@testable import FMIndex
@testable import Alignment
import Htslib

// MARK: - Helpers

/// A parsed SAM alignment record for comparison.
struct ParsedSAMRecord: Sendable {
    let qname: String
    let flag: Int
    let rname: String
    let pos: Int
    let mapq: Int
    let cigar: String
    let rnext: String
    let pnext: Int
    let tlen: Int
    let seq: String
    let tags: [String: String]  // tag -> value

    var isUnmapped: Bool { flag & 0x4 != 0 }
    var isReverse: Bool { flag & 0x10 != 0 }
    var isPaired: Bool { flag & 0x1 != 0 }
    var isProperPair: Bool { flag & 0x2 != 0 }
    var isMateUnmapped: Bool { flag & 0x8 != 0 }
    var isRead1: Bool { flag & 0x40 != 0 }
    var isRead2: Bool { flag & 0x80 != 0 }

    var nm: Int? {
        guard let val = tags["NM"] else { return nil }
        // NM:i:3 or just "3"
        if val.hasPrefix("i:") {
            return Int(val.dropFirst(2))
        }
        return Int(val)
    }

    var md: String? {
        guard let val = tags["MD"] else { return nil }
        if val.hasPrefix("Z:") {
            return String(val.dropFirst(2))
        }
        return val
    }
}

/// Parse a SAM alignment line into a record.
func parseSAMLine(_ line: String) -> ParsedSAMRecord? {
    let fields = line.split(separator: "\t", omittingEmptySubsequences: false)
    guard fields.count >= 11 else { return nil }

    var tags: [String: String] = [:]
    for i in 11..<fields.count {
        let parts = fields[i].split(separator: ":", maxSplits: 1)
        if parts.count == 2 {
            tags[String(parts[0])] = String(parts[1])
        }
    }

    return ParsedSAMRecord(
        qname: String(fields[0]),
        flag: Int(fields[1]) ?? 0,
        rname: String(fields[2]),
        pos: Int(fields[3]) ?? 0,
        mapq: Int(fields[4]) ?? 0,
        cigar: String(fields[5]),
        rnext: String(fields[6]),
        pnext: Int(fields[7]) ?? 0,
        tlen: Int(fields[8]) ?? 0,
        seq: String(fields[9]),
        tags: tags
    )
}

/// Parse a SAM file (alignment records only, no headers) into a dictionary keyed by QNAME.
/// For paired-end, keys include "/1" and "/2" suffixes.
func parseExpectedSAM(path: String) throws -> [String: ParsedSAMRecord] {
    let content = try String(contentsOfFile: path, encoding: .utf8)
    var records: [String: ParsedSAMRecord] = [:]
    for line in content.split(separator: "\n") {
        let lineStr = String(line)
        guard !lineStr.hasPrefix("@"), let record = parseSAMLine(lineStr) else { continue }
        var key = record.qname
        if record.isPaired {
            key += record.isRead1 ? "/1" : "/2"
        }
        records[key] = record
    }
    return records
}

/// Read a FASTQ file and return ReadSequence array.
func loadFASTQ(path: String) throws -> [ReadSequence] {
    let content = try String(contentsOfFile: path, encoding: .utf8)
    let lines = content.split(separator: "\n", omittingEmptySubsequences: false)
    var reads: [ReadSequence] = []
    var i = 0
    while i + 3 < lines.count {
        guard lines[i].hasPrefix("@") else { i += 1; continue }
        let name = String(lines[i].dropFirst())  // strip @
        let seq = String(lines[i + 1])
        let qual = String(lines[i + 3])
        reads.append(ReadSequence(name: name, sequence: seq, qualityString: qual))
        i += 4
    }
    return reads
}

/// Run swiftbwa single-end alignment pipeline and return output SAM records.
func runSingleEndPipeline(
    indexPrefix: String,
    fastqPath: String
) async throws -> [String: ParsedSAMRecord] {
    let index = try FMIndexLoader.load(from: indexPrefix)
    let reads = try loadFASTQ(path: fastqPath)
    let aligner = BWAMemAligner(index: index)

    let tmpPath = NSTemporaryDirectory() + "swiftbwa_e2e_se_\(UUID().uuidString).sam"
    defer { try? FileManager.default.removeItem(atPath: tmpPath) }

    let header = try SAMOutputBuilder.buildHeader(metadata: index.metadata)
    let outFile = try HTSFile(path: tmpPath, mode: "w")
    try header.write(to: outFile)

    try await aligner.alignBatch(reads: reads, outputFile: outFile, header: header)

    // Force flush by dropping the file handle
    _ = consume outFile

    // Read back and parse
    let content = try String(contentsOfFile: tmpPath, encoding: .utf8)
    var records: [String: ParsedSAMRecord] = [:]
    for line in content.split(separator: "\n") {
        let lineStr = String(line)
        guard !lineStr.hasPrefix("@"), let record = parseSAMLine(lineStr) else { continue }
        records[record.qname] = record
    }
    return records
}

/// Run swiftbwa paired-end alignment pipeline and return output SAM records.
func runPairedEndPipeline(
    indexPrefix: String,
    r1Path: String,
    r2Path: String
) async throws -> [String: ParsedSAMRecord] {
    let index = try FMIndexLoader.load(from: indexPrefix)
    let reads1 = try loadFASTQ(path: r1Path)
    let reads2 = try loadFASTQ(path: r2Path)
    var options = BWAMemOptions()
    options.isPairedEnd = true
    let aligner = BWAMemAligner(index: index, options: options)

    let tmpPath = NSTemporaryDirectory() + "swiftbwa_e2e_pe_\(UUID().uuidString).sam"
    defer { try? FileManager.default.removeItem(atPath: tmpPath) }

    let header = try SAMOutputBuilder.buildHeader(metadata: index.metadata)
    let outFile = try HTSFile(path: tmpPath, mode: "w")
    try header.write(to: outFile)

    try await aligner.alignPairedBatch(
        reads1: reads1, reads2: reads2,
        outputFile: outFile, header: header
    )

    _ = consume outFile

    // Read back and parse; key by qname + /1 or /2
    let content = try String(contentsOfFile: tmpPath, encoding: .utf8)
    var records: [String: ParsedSAMRecord] = [:]
    for line in content.split(separator: "\n") {
        let lineStr = String(line)
        guard !lineStr.hasPrefix("@"), let record = parseSAMLine(lineStr) else { continue }
        var key = record.qname
        if record.isPaired {
            key += record.isRead1 ? "/1" : "/2"
        }
        records[key] = record
    }
    return records
}

/// Compare essential FLAG bits between actual and expected.
func flagBitsMatch(_ actual: Int, _ expected: Int) -> Bool {
    let mask = 0x4 | 0x10 | 0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80 | 0x100 | 0x800
    return (actual & mask) == (expected & mask)
}

// MARK: - Test Data Paths

private let testDataDir: String = {
    let file = URL(fileURLWithPath: #filePath)
    let projectRoot = file
        .deletingLastPathComponent()  // BWAMemTests/
        .deletingLastPathComponent()  // Tests/
        .deletingLastPathComponent()  // swiftbwa/
    return projectRoot.appendingPathComponent("TestData").path
}()

private let indexPrefix = testDataDir + "/lambda.fa"

private let indexAvailable = FileManager.default.fileExists(
    atPath: indexPrefix + ".bwt.2bit.64"
)

// MARK: - E2E Single-End Tests

@Suite("E2E Single-End Tests")
struct E2ESingleEndTests {

    @Test("Perfect forward-strand match")
    func testPerfectForwardMatch() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runSingleEndPipeline(
            indexPrefix: indexPrefix,
            fastqPath: testDataDir + "/e2e_reads_se.fq"
        )

        let rec = try #require(results["perfect_fwd_1000"])
        #expect(rec.pos == 1001)
        #expect(rec.cigar == "100M")
        #expect(rec.nm == 0)
        #expect(rec.mapq == 60)
        #expect(!rec.isReverse)
        #expect(!rec.isUnmapped)
    }

    @Test("Second perfect forward-strand match")
    func testPerfectForwardMatch2() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runSingleEndPipeline(
            indexPrefix: indexPrefix,
            fastqPath: testDataDir + "/e2e_reads_se.fq"
        )

        let rec = try #require(results["perfect_fwd_5000"])
        #expect(rec.pos == 5001)
        #expect(rec.cigar == "100M")
        #expect(rec.nm == 0)
        #expect(rec.mapq == 60)
        #expect(!rec.isReverse)
    }

    @Test("Perfect reverse-strand match")
    func testPerfectReverseMatch() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runSingleEndPipeline(
            indexPrefix: indexPrefix,
            fastqPath: testDataDir + "/e2e_reads_se.fq"
        )

        let rec = try #require(results["perfect_rev_10000"])
        #expect(rec.pos == 10001)
        #expect(rec.cigar == "100M")
        #expect(rec.nm == 0)
        #expect(rec.mapq == 60)
        #expect(rec.isReverse)
    }

    @Test("Read with 2 mismatches")
    func testReadWithMismatches() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runSingleEndPipeline(
            indexPrefix: indexPrefix,
            fastqPath: testDataDir + "/e2e_reads_se.fq"
        )

        let rec = try #require(results["snp_2x_20000"])
        #expect(rec.pos == 20001)
        #expect(rec.cigar == "100M")
        #expect(rec.nm == 2)
        #expect(!rec.isReverse)
        // MD tag should contain mismatch markers
        let md = try #require(rec.md)
        #expect(md.contains("A") || md.contains("C") || md.contains("G") || md.contains("T"))
    }

    @Test("Read with 3bp deletion")
    func testReadWithDeletion() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runSingleEndPipeline(
            indexPrefix: indexPrefix,
            fastqPath: testDataDir + "/e2e_reads_se.fq"
        )

        let rec = try #require(results["del3_30000"])
        #expect(rec.pos == 30001)
        // CIGAR should contain a deletion operation
        #expect(rec.cigar.contains("D"))
        #expect((rec.nm ?? 0) >= 3)
        #expect(!rec.isReverse)
    }

    @Test("Unmapped read")
    func testUnmappedRead() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runSingleEndPipeline(
            indexPrefix: indexPrefix,
            fastqPath: testDataDir + "/e2e_reads_se.fq"
        )

        let rec = try #require(results["unmappable"])
        #expect(rec.isUnmapped)
        #expect(rec.pos == 0)
        #expect(rec.cigar == "*")
    }

    @Test("Single-end results match bwa-mem2")
    func testSingleEndVsBwaMem2() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let actual = try await runSingleEndPipeline(
            indexPrefix: indexPrefix,
            fastqPath: testDataDir + "/e2e_reads_se.fq"
        )
        let expected = try parseExpectedSAM(path: testDataDir + "/e2e_expected_se.sam")

        for (name, exp) in expected {
            let act = try #require(actual[name], "Missing read: \(name)")

            // POS: exact match
            #expect(act.pos == exp.pos, "POS mismatch for \(name): \(act.pos) vs \(exp.pos)")

            // FLAG: essential bits match
            #expect(
                flagBitsMatch(act.flag, exp.flag),
                "FLAG mismatch for \(name): \(act.flag) vs \(exp.flag)"
            )

            // CIGAR: exact match for mapped reads
            if !exp.isUnmapped {
                #expect(act.cigar == exp.cigar, "CIGAR mismatch for \(name): \(act.cigar) vs \(exp.cigar)")
            }

            // MAPQ: exact for unique reads (MAPQ=60), ±5 tolerance otherwise
            if exp.mapq == 60 {
                #expect(act.mapq == 60, "MAPQ mismatch for unique read \(name): \(act.mapq)")
            } else {
                #expect(
                    abs(act.mapq - exp.mapq) <= 5,
                    "MAPQ mismatch for \(name): \(act.mapq) vs \(exp.mapq)"
                )
            }

            // NM tag: exact match
            if let expNM = exp.nm {
                #expect(act.nm == expNM, "NM mismatch for \(name): \(String(describing: act.nm)) vs \(expNM)")
            }

            // MD tag: exact match
            if let expMD = exp.md {
                #expect(act.md == expMD, "MD mismatch for \(name): \(String(describing: act.md)) vs \(expMD)")
            }
        }
    }
}

// MARK: - E2E Paired-End Tests

@Suite("E2E Paired-End Tests")
struct E2EPairedEndTests {

    @Test("Proper pair has correct flags and TLEN")
    func testProperPair() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runPairedEndPipeline(
            indexPrefix: indexPrefix,
            r1Path: testDataDir + "/e2e_reads_r1.fq",
            r2Path: testDataDir + "/e2e_reads_r2.fq"
        )

        // Check first proper pair
        let r1 = try #require(results["proper_pair_2000/1"])
        let r2 = try #require(results["proper_pair_2000/2"])

        // Both should be paired
        #expect(r1.isPaired)
        #expect(r2.isPaired)

        // R1 should be read1, R2 should be read2
        #expect(r1.isRead1)
        #expect(r2.isRead2)

        // Neither should be unmapped
        #expect(!r1.isUnmapped)
        #expect(!r2.isUnmapped)

        // Correct positions
        #expect(r1.pos == 2001)
        #expect(r2.pos == 2201)

        // TLEN should be consistent (absolute values match, signs opposite)
        #expect(r1.tlen == -r2.tlen)
        #expect(abs(r1.tlen) == 300)
    }

    @Test("One mate unmapped has correct flags")
    func testOneMateUnmapped() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let results = try await runPairedEndPipeline(
            indexPrefix: indexPrefix,
            r1Path: testDataDir + "/e2e_reads_r1.fq",
            r2Path: testDataDir + "/e2e_reads_r2.fq"
        )

        let r1 = try #require(results["one_unmapped_25000/1"])
        let r2 = try #require(results["one_unmapped_25000/2"])

        // R1 is mapped, R2 is unmapped
        #expect(!r1.isUnmapped)
        #expect(r2.isUnmapped)

        // R1 should have mate-unmapped flag
        #expect(r1.isMateUnmapped)

        // R1 should be at correct position
        #expect(r1.pos == 25001)
        #expect(r1.cigar == "100M")
    }

    @Test("Paired-end results match bwa-mem2")
    func testPairedEndVsBwaMem2() async throws {
        try #require(indexAvailable, "Lambda index not available")

        let actual = try await runPairedEndPipeline(
            indexPrefix: indexPrefix,
            r1Path: testDataDir + "/e2e_reads_r1.fq",
            r2Path: testDataDir + "/e2e_reads_r2.fq"
        )
        let expected = try parseExpectedSAM(path: testDataDir + "/e2e_expected_pe.sam")

        for (name, exp) in expected {
            let act = try #require(actual[name], "Missing read: \(name)")

            // POS: exact match for mapped reads; skip for unmapped (SAM convention difference)
            if !exp.isUnmapped {
                #expect(act.pos == exp.pos, "POS mismatch for \(name): \(act.pos) vs \(exp.pos)")
            }

            // FLAG: essential bits
            #expect(
                flagBitsMatch(act.flag, exp.flag),
                "FLAG mismatch for \(name): \(act.flag) vs \(exp.flag)"
            )

            // CIGAR: exact for mapped reads
            if !exp.isUnmapped {
                #expect(act.cigar == exp.cigar, "CIGAR mismatch for \(name): \(act.cigar) vs \(exp.cigar)")
            }

            // MAPQ: exact for unique (60), ±5 otherwise
            if exp.mapq == 60 {
                #expect(act.mapq == 60, "MAPQ mismatch for unique \(name): \(act.mapq)")
            } else if exp.mapq > 0 {
                #expect(
                    abs(act.mapq - exp.mapq) <= 5,
                    "MAPQ mismatch for \(name): \(act.mapq) vs \(exp.mapq)"
                )
            }

            // NM tag: exact match
            if let expNM = exp.nm {
                #expect(act.nm == expNM, "NM mismatch for \(name): \(String(describing: act.nm)) vs \(expNM)")
            }

            // MD tag: exact match
            if let expMD = exp.md {
                #expect(act.md == expMD, "MD mismatch for \(name): \(String(describing: act.md)) vs \(expMD)")
            }
        }
    }
}
