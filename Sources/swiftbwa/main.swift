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

    // Algorithm options
    @Option(name: .shortAndLong, help: "Number of threads")
    var threads: Int = 1

    @Option(name: .short, help: "Minimum seed length")
    var k: Int32 = 19

    @Option(name: .short, help: "Band width")
    var w: Int32 = 100

    @Option(name: .short, help: "Off-diagonal X-dropoff")
    var d: Int32 = 100

    @Option(name: .short, help: "Seed split ratio")
    var r: Float = 1.5

    @Option(name: .short, help: "Skip seeds with more than INT occurrences")
    var c: Int32 = 500

    @Option(name: [.customShort("D")], help: "Drop chains shorter than FLOAT fraction of longest overlapping chain")
    var chainDrop: Float = 0.50

    @Option(name: [.customShort("W")], help: "Discard chains with seeded bases shorter than INT")
    var minChainWeight: Int32 = 0

    @Option(name: .short, help: "Max mate rescue rounds per read")
    var m: Int32 = 50

    @Flag(name: [.customShort("S")], help: "Skip mate rescue")
    var skipRescue: Bool = false

    @Flag(name: [.customShort("P")], help: "Skip pairing; mate rescue still performed unless -S")
    var skipPairing: Bool = false

    // Scoring options
    @Option(name: [.customShort("A")], help: "Score for a sequence match (scales -TdBOELU)")
    var matchScore: Int32 = 1

    @Option(name: [.customShort("B")], help: "Mismatch penalty")
    var mismatch: Int32 = 4

    @Option(name: [.customShort("O")], help: "Gap open penalty [,INT for deletions]")
    var gapOpen: String = "6"

    @Option(name: [.customShort("E")], help: "Gap extension penalty [,INT for deletions]")
    var gapExtend: String = "1"

    @Option(name: [.customShort("L")], help: "Clipping penalty [,INT for 3' end]")
    var clipPenalty: String = "5"

    @Option(name: [.customShort("U")], help: "Penalty for unpaired read pair")
    var unpairedPenalty: Int32 = 17

    @Option(name: [.customShort("T")], help: "Minimum alignment score to output")
    var minScore: Int32 = 30

    @Option(name: [.customLong("xa-hits")], help: "XA hit limits [,INT for ALT contigs]")
    var xaHits: String = "5,200"

    // Input/output options
    @Option(name: [.customShort("R")], help: "Read group header line (@RG\\tID:...)")
    var readGroup: String?

    @Option(name: [.customShort("H")], help: "Insert header lines (string starting with @ or file path)")
    var headerLines: String?

    @Option(name: .shortAndLong, help: "Output file [stdout]")
    var output: String?

    @Flag(name: [.customShort("M")], help: "Mark shorter split hits as secondary")
    var markSecondary: Bool = false

    @Flag(name: [.customShort("Y")], help: "Use soft clipping for supplementary alignments")
    var softClipSupp: Bool = false

    @Flag(name: [.customShort("5")], help: "Take split alignment with smallest coordinate as primary")
    var primary5: Bool = false

    @Flag(name: .short, help: "Don't modify MAPQ of supplementary alignments")
    var q: Bool = false

    @Flag(name: .short, help: "Treat ALT contigs as part of primary assembly")
    var j: Bool = false

    @Flag(name: [.customShort("a")], help: "Output all alignments for SE or unpaired PE")
    var outputAll: Bool = false

    @Argument(help: "Index prefix (reference genome)")
    var indexPrefix: String

    @Argument(help: "FASTQ file(s)")
    var fastqFiles: [String]

    func run() async throws {
        // Parse comma-separated penalty pairs
        let gapOpenParts = gapOpen.split(separator: ",").compactMap { Int32($0) }
        let gapExtendParts = gapExtend.split(separator: ",").compactMap { Int32($0) }
        let clipParts = clipPenalty.split(separator: ",").compactMap { Int32($0) }
        let xaParts = xaHits.split(separator: ",").compactMap { Int32($0) }

        // Build scoring parameters
        var scoring = ScoringParameters()
        scoring.matchScore = matchScore
        scoring.minSeedLength = k
        scoring.bandWidth = w
        scoring.zDrop = d
        scoring.seedSplitRatio = r
        scoring.maxOccurrences = c
        scoring.chainDropRatio = chainDrop
        scoring.minChainWeight = minChainWeight
        scoring.maxMatesw = m
        scoring.unpairedPenalty = unpairedPenalty
        scoring.numThreads = threads

        // Scale penalties by match score if -A != 1 and user didn't override
        let a = matchScore
        scoring.mismatchPenalty = mismatch * a
        scoring.gapOpenPenalty = gapOpenParts[0] * a
        scoring.gapExtendPenalty = gapExtendParts[0] * a
        scoring.gapOpenPenaltyDeletion = (gapOpenParts.count > 1 ? gapOpenParts[1] : gapOpenParts[0]) * a
        scoring.gapExtendPenaltyDeletion = (gapExtendParts.count > 1 ? gapExtendParts[1] : gapExtendParts[0]) * a
        scoring.penClip5 = (clipParts.first ?? 5) * a
        scoring.penClip3 = (clipParts.count > 1 ? clipParts[1] : clipParts.first ?? 5) * a
        scoring.minOutputScore = minScore * a
        scoring.zDrop = d * a

        scoring.maxXAHits = xaParts.first ?? 5
        scoring.maxXAHitsAlt = xaParts.count > 1 ? xaParts[1] : 200

        // Set behavioral flags
        var flagBits: Int32 = 0
        if markSecondary { flagBits |= ScoringParameters.flagNoMulti }
        if softClipSupp { flagBits |= ScoringParameters.flagSoftClip }
        if primary5 { flagBits |= ScoringParameters.flagPrimary5 | ScoringParameters.flagKeepSuppMapq }
        if q { flagBits |= ScoringParameters.flagKeepSuppMapq }
        if skipRescue { flagBits |= ScoringParameters.flagNoRescue }
        if skipPairing { flagBits |= ScoringParameters.flagNoPairing }
        if j { flagBits |= ScoringParameters.flagNoAlt }
        if outputAll { flagBits |= ScoringParameters.flagAll }
        scoring.flag = flagBits

        var options = BWAMemOptions()
        options.scoring = scoring
        options.isPairedEnd = fastqFiles.count >= 2
        options.readGroupLine = readGroup
        options.ignoreAlt = j

        // Handle -H: string starting with @ or file path
        if let hArg = headerLines {
            if hArg.hasPrefix("@") {
                options.headerLines = hArg
            } else {
                options.headerLines = try String(contentsOfFile: hArg, encoding: .utf8)
            }
        }

        // Load index
        fputs("Loading index from \(indexPrefix)...\n", stderr)
        let index = try FMIndexLoader.load(from: indexPrefix, skipAlt: options.ignoreAlt)
        fputs("Index loaded: ref_len=\(index.referenceSeqLen), "
              + "genome_len=\(index.genomeLength), "
              + "\(index.metadata.numSequences) sequences\n", stderr)

        // Output
        let outputMode = output != nil ? "wb" : "w"
        let outputPath = output ?? "-"

        let aligner = BWAMemAligner(index: index, options: options)
        let header = try SAMOutputBuilder.buildHeader(
            metadata: index.metadata,
            readGroupLine: options.readGroupLine,
            headerLines: options.headerLines
        )
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
