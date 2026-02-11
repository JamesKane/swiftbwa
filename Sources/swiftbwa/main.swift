import ArgumentParser
import BWAMem
import FMIndex
import BWACore
import Foundation
import Htslib

@main
struct SwiftBWA: AsyncParsableCommand {
    static let configuration = CommandConfiguration(
        commandName: "swiftbwa",
        abstract: "BWA-MEM sequence aligner reimplemented in Swift",
        version: "0.1.0",
        subcommands: [Mem.self],
        defaultSubcommand: Mem.self
    )
}

struct Mem: AsyncParsableCommand {
    static let configuration = CommandConfiguration(
        abstract: "Align reads using BWA-MEM algorithm"
    )

    @Option(name: .shortAndLong, help: "Number of threads")
    var threads: Int = 1

    @Option(name: .short, help: "Minimum seed length")
    var k: Int32 = 19

    @Option(name: .short, help: "Band width")
    var w: Int32 = 100

    @Option(name: .short, help: "Off-diagonal X-dropoff")
    var d: Int32 = 100

    @Option(name: [.customShort("B")], help: "Mismatch penalty")
    var mismatch: Int32 = 4

    @Option(name: [.customShort("O")], help: "Gap open penalty")
    var gapOpen: Int32 = 6

    @Option(name: [.customShort("E")], help: "Gap extension penalty")
    var gapExtend: Int32 = 1

    @Option(name: [.customShort("T")], help: "Minimum alignment score to output")
    var minScore: Int32 = 30

    @Option(name: .shortAndLong, help: "Output file [stdout]")
    var output: String?

    @Argument(help: "Index prefix (reference genome)")
    var indexPrefix: String

    @Argument(help: "FASTQ file(s)")
    var fastqFiles: [String]

    func run() async throws {
        // Build scoring parameters
        var scoring = ScoringParameters()
        scoring.minSeedLength = k
        scoring.bandWidth = w
        scoring.zDrop = d
        scoring.mismatchPenalty = mismatch
        scoring.gapOpenPenalty = gapOpen
        scoring.gapExtendPenalty = gapExtend
        scoring.gapOpenPenaltyDeletion = gapOpen
        scoring.gapExtendPenaltyDeletion = gapExtend
        scoring.minOutputScore = minScore
        scoring.numThreads = threads

        var options = BWAMemOptions()
        options.scoring = scoring
        options.isPairedEnd = fastqFiles.count >= 2

        // Load index
        fputs("Loading index from \(indexPrefix)...\n", stderr)
        let index = try FMIndexLoader.load(from: indexPrefix)
        fputs("Index loaded: ref_len=\(index.referenceSeqLen), "
              + "genome_len=\(index.genomeLength), "
              + "\(index.metadata.numSequences) sequences\n", stderr)

        // Parse FASTQ and align
        let reads = try parseFASTQ(path: fastqFiles[0])
        fputs("Loaded \(reads.count) reads from \(fastqFiles[0])\n", stderr)

        // Output
        let outputMode = output != nil ? "wb" : "w"
        let outputPath = output ?? "-"

        fputs("Aligning...\n", stderr)

        let aligner = BWAMemAligner(index: index, options: options)
        let header = try SAMOutputBuilder.buildHeader(metadata: index.metadata)
        let outFile = try HTSFile(path: outputPath, mode: outputMode)
        try header.write(to: outFile)

        try await aligner.alignBatch(
            reads: reads,
            outputFile: outFile,
            header: header
        )

        fputs("Done.\n", stderr)
    }
}

// MARK: - Simple FASTQ Parser

func parseFASTQ(path: String) throws -> [ReadSequence] {
    guard let data = FileManager.default.contents(atPath: path),
          let content = String(data: data, encoding: .utf8) else {
        throw BWAError.fileIOError("Cannot read FASTQ file: \(path)")
    }

    var reads: [ReadSequence] = []
    let lines = content.split(separator: "\n", omittingEmptySubsequences: false)
    var i = 0

    while i + 3 < lines.count {
        let header = lines[i]
        let sequence = lines[i + 1]
        let qualStr = lines[i + 3]

        guard header.hasPrefix("@") else {
            i += 1
            continue
        }

        let name = String(header.dropFirst().split(separator: " ", maxSplits: 1).first ?? "")
        let comment = header.contains(" ")
            ? String(header.dropFirst().split(separator: " ", maxSplits: 1).last ?? "")
            : ""

        let read = ReadSequence(
            name: name,
            sequence: String(sequence),
            qualityString: String(qualStr),
            comment: comment
        )
        reads.append(read)
        i += 4
    }

    return reads
}
