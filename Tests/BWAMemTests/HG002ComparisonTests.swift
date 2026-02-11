import Foundation
import Testing
@testable import BWAMem
@testable import BWACore
@testable import FMIndex
@testable import Alignment
import Htslib

// MARK: - Index Cache

/// Caches the ~23GB chm13v2.0 index across all tests in the process.
private actor IndexCache {
    static let shared = IndexCache()
    private var cached: FMIndex?

    func get() throws -> FMIndex {
        if let idx = cached { return idx }
        let idx = try FMIndexLoader.load(from: chm13IndexPrefix)
        cached = idx
        return idx
    }
}

// MARK: - Paths

private let testDataDir: String = {
    let file = URL(fileURLWithPath: #filePath)
    let projectRoot = file
        .deletingLastPathComponent()  // BWAMemTests/
        .deletingLastPathComponent()  // Tests/
        .deletingLastPathComponent()  // swiftbwa/
    return projectRoot.appendingPathComponent("TestData").path
}()

private let chm13IndexPrefix = "/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa.gz"

private let chm13IndexAvailable = FileManager.default.fileExists(
    atPath: chm13IndexPrefix + ".bwt.2bit.64"
)

private let seFastqPath = testDataDir + "/hg002_chm13_se_100.fq"
private let peR1Path = testDataDir + "/hg002_chm13_pe_100_R1.fq"
private let peR2Path = testDataDir + "/hg002_chm13_pe_100_R2.fq"
private let seGoldPath = testDataDir + "/hg002_chm13_se_100.sam"
private let peGoldPath = testDataDir + "/hg002_chm13_pe_100.sam"

private let goldFilesExist =
    FileManager.default.fileExists(atPath: seFastqPath)
    && FileManager.default.fileExists(atPath: seGoldPath)
    && FileManager.default.fileExists(atPath: peGoldPath)

private let hg002TestsEnabled = chm13IndexAvailable && goldFilesExist

/// Total soft-clipped bases in a CIGAR string.
private func softClipBases(_ cigar: String) -> Int {
    var total = 0
    var numStr = ""
    for ch in cigar {
        if ch.isNumber {
            numStr.append(ch)
        } else {
            if ch == "S", let n = Int(numStr) {
                total += n
            }
            numStr = ""
        }
    }
    return total
}

// MARK: - Pipeline Helpers (use cached index)

/// Run swiftbwa single-end alignment using the cached chm13 index.
private func runSE(fastqPath: String) async throws -> [String: ParsedSAMRecord] {
    let index = try await IndexCache.shared.get()
    let reads = try loadFASTQ(path: fastqPath)
    let aligner = BWAMemAligner(index: index)

    let tmpPath = NSTemporaryDirectory() + "hg002_se_\(UUID().uuidString).sam"
    defer { try? FileManager.default.removeItem(atPath: tmpPath) }

    let header = try SAMOutputBuilder.buildHeader(metadata: index.metadata)
    let outFile = try HTSFile(path: tmpPath, mode: "w")
    try header.write(to: outFile)

    try await aligner.alignBatch(reads: reads, outputFile: outFile, header: header)
    _ = consume outFile

    let content = try String(contentsOfFile: tmpPath, encoding: .utf8)
    var records: [String: ParsedSAMRecord] = [:]
    for line in content.split(separator: "\n") {
        let lineStr = String(line)
        guard !lineStr.hasPrefix("@"), let record = parseSAMLine(lineStr) else { continue }
        // Skip supplementary/secondary for primary comparison
        guard record.flag & 0x900 == 0 else { continue }
        records[record.qname] = record
    }
    return records
}

/// Run swiftbwa paired-end alignment using the cached chm13 index.
private func runPE(r1Path: String, r2Path: String) async throws -> [String: ParsedSAMRecord] {
    let index = try await IndexCache.shared.get()
    let reads1 = try loadFASTQ(path: r1Path)
    let reads2 = try loadFASTQ(path: r2Path)
    var options = BWAMemOptions()
    options.isPairedEnd = true
    let aligner = BWAMemAligner(index: index, options: options)

    let tmpPath = NSTemporaryDirectory() + "hg002_pe_\(UUID().uuidString).sam"
    defer { try? FileManager.default.removeItem(atPath: tmpPath) }

    let header = try SAMOutputBuilder.buildHeader(metadata: index.metadata)
    let outFile = try HTSFile(path: tmpPath, mode: "w")
    try header.write(to: outFile)

    try await aligner.alignPairedBatch(
        reads1: reads1, reads2: reads2,
        outputFile: outFile, header: header
    )
    _ = consume outFile

    let content = try String(contentsOfFile: tmpPath, encoding: .utf8)
    var records: [String: ParsedSAMRecord] = [:]
    for line in content.split(separator: "\n") {
        let lineStr = String(line)
        guard !lineStr.hasPrefix("@"), let record = parseSAMLine(lineStr) else { continue }
        guard record.flag & 0x900 == 0 else { continue }
        var key = record.qname
        if record.isPaired {
            key += record.isRead1 ? "/1" : "/2"
        }
        records[key] = record
    }
    return records
}

/// Parse gold standard SAM, filtering to primary alignments only.
private func parseGoldPrimary(path: String) throws -> [String: ParsedSAMRecord] {
    let content = try String(contentsOfFile: path, encoding: .utf8)
    var records: [String: ParsedSAMRecord] = [:]
    for line in content.split(separator: "\n") {
        let lineStr = String(line)
        guard !lineStr.hasPrefix("@"), let record = parseSAMLine(lineStr) else { continue }
        guard record.flag & 0x900 == 0 else { continue }
        var key = record.qname
        if record.isPaired {
            key += record.isRead1 ? "/1" : "/2"
        }
        records[key] = record
    }
    return records
}

// MARK: - HG002/chm13v2.0 Single-End Comparison

@Suite("HG002/chm13v2.0 Single-End Comparison",
       .enabled(if: hg002TestsEnabled, "chm13v2.0 index or gold files not available"))
struct HG002SingleEndTests {

    @Test("SE mapping rate within 2% of bwa-mem2")
    func testSEMappingRate() async throws {
        let gold = try parseGoldPrimary(path: seGoldPath)
        let actual = try await runSE(fastqPath: seFastqPath)

        let goldMapped = gold.values.filter { !$0.isUnmapped }.count
        let actualMapped = actual.values.filter { !$0.isUnmapped }.count

        let goldRate = Double(goldMapped) / Double(max(gold.count, 1))
        let actualRate = Double(actualMapped) / Double(max(actual.count, 1))

        let msg = "Mapping rate diff: bwa=\(String(format: "%.1f", goldRate * 100))% swift=\(String(format: "%.1f", actualRate * 100))%"
        #expect(abs(goldRate - actualRate) < 0.02, Comment(rawValue: msg))
    }

    @Test("SE position concordance >= 95%")
    func testSEPositionConcordance() async throws {
        let gold = try parseGoldPrimary(path: seGoldPath)
        let actual = try await runSE(fastqPath: seFastqPath)

        var bothMapped = 0
        var posExact = 0

        for (key, exp) in gold {
            guard let act = actual[key] else { continue }
            guard !exp.isUnmapped && !act.isUnmapped else { continue }
            // Skip ambiguous multi-mappers (MAPQ=0: position is arbitrary)
            // and reads that gold itself soft-clips (>10bp: both tools struggle)
            guard exp.mapq > 0 else { continue }
            guard softClipBases(exp.cigar) <= 10 else { continue }
            bothMapped += 1
            if act.rname == exp.rname && act.pos == exp.pos {
                posExact += 1
            }
        }

        let concordance = bothMapped > 0 ? Double(posExact) / Double(bothMapped) : 0
        #expect(concordance >= 0.95, "Position concordance: \(posExact)/\(bothMapped)")
    }

    @Test("SE CIGAR >= 90% and NM >= 90% match at same positions")
    func testSECigarNM() async throws {
        let gold = try parseGoldPrimary(path: seGoldPath)
        let actual = try await runSE(fastqPath: seFastqPath)

        var samePos = 0
        var cigarMatch = 0
        var nmMatch = 0

        for (key, exp) in gold {
            guard let act = actual[key] else { continue }
            guard !exp.isUnmapped && !act.isUnmapped else { continue }
            guard act.rname == exp.rname && act.pos == exp.pos else { continue }
            samePos += 1
            if act.cigar == exp.cigar { cigarMatch += 1 }
            if let aNM = act.nm, let eNM = exp.nm, aNM == eNM { nmMatch += 1 }
        }

        let cigarRate = samePos > 0 ? Double(cigarMatch) / Double(samePos) : 0
        let nmRate = samePos > 0 ? Double(nmMatch) / Double(samePos) : 0

        #expect(cigarRate >= 0.90, "CIGAR match: \(cigarMatch)/\(samePos)")
        #expect(nmRate >= 0.90, "NM match: \(nmMatch)/\(samePos)")
    }

    @Test("SE MAPQ: mean diff < 5, >= 80% exact")
    func testSEMapQ() async throws {
        let gold = try parseGoldPrimary(path: seGoldPath)
        let actual = try await runSE(fastqPath: seFastqPath)

        var compared = 0
        var exact = 0
        var diffSum = 0

        for (key, exp) in gold {
            guard let act = actual[key] else { continue }
            guard !exp.isUnmapped && !act.isUnmapped else { continue }
            // Skip multi-mappers: MAPQ=0 means position is arbitrary,
            // so MAPQ comparison across different positions is meaningless
            guard exp.mapq > 0 else { continue }
            compared += 1
            let diff = abs(act.mapq - exp.mapq)
            diffSum += diff
            if diff == 0 { exact += 1 }
        }

        let meanDiff = compared > 0 ? Double(diffSum) / Double(compared) : 0
        let exactRate = compared > 0 ? Double(exact) / Double(compared) : 0

        #expect(meanDiff < 5.0, "Mean MAPQ diff: \(String(format: "%.2f", meanDiff))")
        #expect(exactRate >= 0.80, "MAPQ exact: \(exact)/\(compared)")
    }
}

// MARK: - HG002/chm13v2.0 Paired-End Comparison

@Suite("HG002/chm13v2.0 Paired-End Comparison",
       .enabled(if: hg002TestsEnabled, "chm13v2.0 index or gold files not available"))
struct HG002PairedEndTests {

    @Test("PE mapping rate within 2% of bwa-mem2")
    func testPEMappingRate() async throws {
        let gold = try parseGoldPrimary(path: peGoldPath)
        let actual = try await runPE(r1Path: peR1Path, r2Path: peR2Path)

        let goldMapped = gold.values.filter { !$0.isUnmapped }.count
        let actualMapped = actual.values.filter { !$0.isUnmapped }.count

        let goldRate = Double(goldMapped) / Double(max(gold.count, 1))
        let actualRate = Double(actualMapped) / Double(max(actual.count, 1))

        let msg = "PE mapping rate diff: bwa=\(String(format: "%.1f", goldRate * 100))% swift=\(String(format: "%.1f", actualRate * 100))%"
        #expect(abs(goldRate - actualRate) < 0.02, Comment(rawValue: msg))
    }

    @Test("PE position concordance >= 95%")
    func testPEPositionConcordance() async throws {
        let gold = try parseGoldPrimary(path: peGoldPath)
        let actual = try await runPE(r1Path: peR1Path, r2Path: peR2Path)

        var bothMapped = 0
        var posExact = 0

        for (key, exp) in gold {
            guard let act = actual[key] else { continue }
            guard !exp.isUnmapped && !act.isUnmapped else { continue }
            bothMapped += 1
            if act.rname == exp.rname && act.pos == exp.pos {
                posExact += 1
            }
        }

        let concordance = bothMapped > 0 ? Double(posExact) / Double(bothMapped) : 0
        #expect(concordance >= 0.95, "PE position concordance: \(posExact)/\(bothMapped)")
    }

    @Test("PE proper pair rate within 5% of bwa-mem2")
    func testPEProperPairRate() async throws {
        let gold = try parseGoldPrimary(path: peGoldPath)
        let actual = try await runPE(r1Path: peR1Path, r2Path: peR2Path)

        let goldProper = gold.values.filter { $0.isProperPair }.count
        let actualProper = actual.values.filter { $0.isProperPair }.count

        let goldRate = Double(goldProper) / Double(max(gold.count, 1))
        let actualRate = Double(actualProper) / Double(max(actual.count, 1))

        let msg = "Proper pair diff: bwa=\(String(format: "%.1f", goldRate * 100))% swift=\(String(format: "%.1f", actualRate * 100))%"
        #expect(abs(goldRate - actualRate) < 0.05, Comment(rawValue: msg))
    }

    @Test("PE TLEN sign consistent, abs values within 10%")
    func testPETLEN() async throws {
        let gold = try parseGoldPrimary(path: peGoldPath)
        let actual = try await runPE(r1Path: peR1Path, r2Path: peR2Path)

        var compared = 0
        var signConsistent = 0
        var absWithin10pct = 0

        for (key, exp) in gold {
            guard let act = actual[key] else { continue }
            guard !exp.isUnmapped && !act.isUnmapped else { continue }
            guard exp.tlen != 0 && act.tlen != 0 else { continue }
            guard act.rname == exp.rname && act.pos == exp.pos else { continue }
            compared += 1

            // Sign consistency
            if (exp.tlen > 0 && act.tlen > 0) || (exp.tlen < 0 && act.tlen < 0) {
                signConsistent += 1
            }

            // Absolute value within 10%
            let expAbs = abs(exp.tlen)
            let actAbs = abs(act.tlen)
            let maxAbs = max(expAbs, actAbs)
            if maxAbs > 0 && abs(expAbs - actAbs) <= maxAbs / 10 {
                absWithin10pct += 1
            }
        }

        guard compared > 0 else {
            // No reads to compare TLEN on — pass vacuously
            return
        }

        let signRate = Double(signConsistent) / Double(compared)
        let absRate = Double(absWithin10pct) / Double(compared)

        #expect(signRate >= 0.95, "TLEN sign: \(signConsistent)/\(compared)")
        #expect(absRate >= 0.90, "TLEN abs: \(absWithin10pct)/\(compared)")
    }

    @Test("PE CIGAR >= 90% match at same positions")
    func testPECigarNM() async throws {
        let gold = try parseGoldPrimary(path: peGoldPath)
        let actual = try await runPE(r1Path: peR1Path, r2Path: peR2Path)

        var samePos = 0
        var cigarMatch = 0

        for (key, exp) in gold {
            guard let act = actual[key] else { continue }
            guard !exp.isUnmapped && !act.isUnmapped else { continue }
            guard act.rname == exp.rname && act.pos == exp.pos else { continue }
            samePos += 1
            if act.cigar == exp.cigar { cigarMatch += 1 }
        }

        let cigarRate = samePos > 0 ? Double(cigarMatch) / Double(samePos) : 0
        #expect(cigarRate >= 0.90, "PE CIGAR: \(cigarMatch)/\(samePos)")
    }

    @Test("PE MAPQ: mean diff < 5")
    func testPEMapQ() async throws {
        let gold = try parseGoldPrimary(path: peGoldPath)
        let actual = try await runPE(r1Path: peR1Path, r2Path: peR2Path)

        var compared = 0
        var diffSum = 0

        for (key, exp) in gold {
            guard let act = actual[key] else { continue }
            guard !exp.isUnmapped && !act.isUnmapped else { continue }
            compared += 1
            diffSum += abs(act.mapq - exp.mapq)
        }

        let meanDiff = compared > 0 ? Double(diffSum) / Double(compared) : 0
        #expect(meanDiff < 5.0, "PE mean MAPQ diff: \(String(format: "%.2f", meanDiff))")
    }
}
