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

        // Output
        let outputMode = output != nil ? "wb" : "w"
        let outputPath = output ?? "-"

        let aligner = BWAMemAligner(index: index, options: options)
        let header = try SAMOutputBuilder.buildHeader(metadata: index.metadata)
        let outFile = try HTSFile(path: outputPath, mode: outputMode)
        try header.write(to: outFile)

        if options.isPairedEnd {
            // Paired-end mode
            let reads1 = try readFASTQ(path: fastqFiles[0])
            let reads2 = try readFASTQ(path: fastqFiles[1])
            guard reads1.count == reads2.count else {
                throw BWAError.fileIOError(
                    "R1 and R2 have different read counts: "
                    + "\(reads1.count) vs \(reads2.count)"
                )
            }
            fputs("Loaded \(reads1.count) paired reads\n", stderr)
            fputs("Aligning paired-end...\n", stderr)

            try await aligner.alignPairedBatch(
                reads1: reads1, reads2: reads2,
                outputFile: outFile, header: header
            )
        } else {
            // Single-end mode
            let reads = try readFASTQ(path: fastqFiles[0])
            fputs("Loaded \(reads.count) reads from \(fastqFiles[0])\n", stderr)
            fputs("Aligning...\n", stderr)

            try await aligner.alignBatch(
                reads: reads,
                outputFile: outFile,
                header: header
            )
        }

        fputs("Done.\n", stderr)
    }
}

// MARK: - FASTQ Reader

func readFASTQ(path: String) throws -> [ReadSequence] {
    let file = try HTSFile(path: path, mode: "r")
    let header = try file.samHeader()
    let iterator = file.samIterator(header: header)

    var reads: [ReadSequence] = []
    while let record = iterator.next() {
        let name = record.queryName
        let seq = record.sequence
        let quals = record.qualities

        let bases: [UInt8] = (0..<seq.count).map { i in
            switch seq[i] {
            case "A", "a": return 0
            case "C", "c": return 1
            case "G", "g": return 2
            case "T", "t": return 3
            default:       return 4  // N
            }
        }
        let qualities: [UInt8] = (0..<quals.count).map { quals[$0] }

        reads.append(ReadSequence(name: name, bases: bases, qualities: qualities))
    }
    return reads
}
