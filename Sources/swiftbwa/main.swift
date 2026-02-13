import ArgumentParser
import BWAMem
import FMIndex
import BWACore
import Foundation
import Htslib
#if canImport(Metal)
import MetalSW
#endif

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

    @Option(name: [.customShort("I")], help: "Insert size distribution: mean,stddev[,max[,min]]")
    var insertSize: String?

    @Option(name: [.customShort("K")], help: "Process INT input bases per batch for reproducibility")
    var batchSize: Int?

    @Option(name: .short, help: "Seed occurrence re-seeding threshold length")
    var y: Int32 = 0

    @Option(name: .short, help: "Verbosity level: 1=error, 2=warning, 3=info, 4+=debug")
    var v: Int = 3

    @Flag(name: [.customShort("C")], help: "Append FASTQ comment to SAM output")
    var appendComment: Bool = false

    @Flag(name: [.customShort("V")], help: "Output reference header in XR tag")
    var outputRefHeader: Bool = false

    @Flag(name: .short, help: "First query file is interleaved paired-end")
    var p: Bool = false

    @Flag(name: .long, help: "Use Metal GPU acceleration for Smith-Waterman")
    var gpu: Bool = false

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

        // Scale penalties by match score (bwa-mem2 only scales when -A != 1)
        let a = matchScore
        if a != 1 {
            scoring.mismatchPenalty = mismatch * a
            scoring.gapOpenPenalty = gapOpenParts[0] * a
            scoring.gapExtendPenalty = gapExtendParts[0] * a
            scoring.gapOpenPenaltyDeletion = (gapOpenParts.count > 1 ? gapOpenParts[1] : gapOpenParts[0]) * a
            scoring.gapExtendPenaltyDeletion = (gapExtendParts.count > 1 ? gapExtendParts[1] : gapExtendParts[0]) * a
            scoring.penClip5 = (clipParts.first ?? 5) * a
            scoring.penClip3 = (clipParts.count > 1 ? clipParts[1] : clipParts.first ?? 5) * a
            scoring.minOutputScore = minScore * a
            scoring.zDrop = d * a
        } else {
            scoring.mismatchPenalty = mismatch
            scoring.gapOpenPenalty = gapOpenParts[0]
            scoring.gapExtendPenalty = gapExtendParts[0]
            scoring.gapOpenPenaltyDeletion = gapOpenParts.count > 1 ? gapOpenParts[1] : gapOpenParts[0]
            scoring.gapExtendPenaltyDeletion = gapExtendParts.count > 1 ? gapExtendParts[1] : gapExtendParts[0]
            scoring.penClip5 = clipParts.first ?? 5
            scoring.penClip3 = clipParts.count > 1 ? clipParts[1] : clipParts.first ?? 5
            scoring.minOutputScore = minScore
            scoring.zDrop = d
        }

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
        options.isPairedEnd = fastqFiles.count >= 2 || p
        options.readGroupLine = readGroup
        options.ignoreAlt = j
        options.appendComment = appendComment
        options.outputRefHeader = outputRefHeader
        options.verbosity = v
        options.useGPU = gpu
        if let bs = batchSize { scoring.chunkSize = bs }
        scoring.reseedLength = y

        // Parse -I insert size override
        if let isStr = insertSize {
            let parts = isStr.split(separator: ",").compactMap { Double($0) }
            guard parts.count >= 2 else {
                throw BWAError.fileIOError("-I requires at least mean,stddev")
            }
            options.manualInsertSize = InsertSizeOverride(
                mean: parts[0],
                stddev: parts[1],
                max: parts.count > 2 ? parts[2] : nil,
                min: parts.count > 3 ? parts[3] : nil
            )
        }

        // Handle -H: string starting with @ or file path
        if let hArg = headerLines {
            if hArg.hasPrefix("@") {
                options.headerLines = hArg
            } else {
                options.headerLines = try String(contentsOfFile: hArg, encoding: .utf8)
            }
        }

        // Load index
        if v >= 3 { fputs("Loading index from \(indexPrefix)...\n", stderr) }
        let index = try FMIndexLoader.load(from: indexPrefix, skipAlt: options.ignoreAlt)
        if v >= 3 {
            fputs("Index loaded: ref_len=\(index.referenceSeqLen), "
                  + "genome_len=\(index.genomeLength), "
                  + "\(index.metadata.numSequences) sequences\n", stderr)
        }

        // Report GPU status
        if gpu && v >= 3 {
            #if canImport(Metal)
            if let engine = MetalSWEngine.shared {
                fputs("[GPU] Metal acceleration enabled: \(engine.device.name)\n", stderr)
            } else {
                fputs("[GPU] Warning: Metal device not available, falling back to CPU\n", stderr)
            }
            #else
            fputs("[GPU] Warning: Metal not available on this platform, falling back to CPU\n", stderr)
            #endif
        }

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

        let chunkBases = batchSize ?? defaultChunkBases

        if p {
            // Interleaved paired-end mode: stream pairs from one file
            let inFile = try HTSFile(path: fastqFiles[0], mode: "r")
            let inHeader = try inFile.samHeader()
            let iter = inFile.samIterator(header: inHeader)
            var totalPairs = 0
            var totalBases = 0
            if v >= 3 { fputs("Aligning interleaved paired-end...\n", stderr) }
            while true {
                let (chunk1, chunk2) = readInterleavedChunk(
                    from: iter, maxBases: chunkBases,
                    preserveComments: appendComment
                )
                if chunk1.isEmpty { break }
                guard chunk1.count == chunk2.count else {
                    throw BWAError.fileIOError(
                        "Interleaved FASTQ has odd number of records"
                    )
                }
                totalPairs += chunk1.count
                totalBases += chunk1.reduce(0) { $0 + $1.length }
                    + chunk2.reduce(0) { $0 + $1.length }
                if v >= 3 {
                    fputs("[M::main] Processed \(totalPairs) pairs "
                        + "(\(String(format: "%.1f", Double(totalBases) / 1_000_000))M bases)\n",
                        stderr)
                }
                try await aligner.alignPairedBatch(
                    reads1: chunk1, reads2: chunk2,
                    outputFile: outFile, header: header
                )
            }
        } else if options.isPairedEnd {
            // Two-file paired-end mode: stream from both files
            let inFile1 = try HTSFile(path: fastqFiles[0], mode: "r")
            let inHeader1 = try inFile1.samHeader()
            let iter1 = inFile1.samIterator(header: inHeader1)
            let inFile2 = try HTSFile(path: fastqFiles[1], mode: "r")
            let inHeader2 = try inFile2.samHeader()
            let iter2 = inFile2.samIterator(header: inHeader2)
            var totalPairs = 0
            var totalBases = 0
            if v >= 3 { fputs("Aligning paired-end...\n", stderr) }
            while true {
                let chunk1 = readChunk(
                    from: iter1, maxBases: chunkBases,
                    preserveComments: appendComment
                )
                if chunk1.isEmpty { break }
                let chunk2 = readExactCount(
                    from: iter2, count: chunk1.count,
                    preserveComments: appendComment
                )
                guard chunk2.count == chunk1.count else {
                    throw BWAError.fileIOError(
                        "R1 and R2 have different read counts"
                    )
                }
                totalPairs += chunk1.count
                totalBases += chunk1.reduce(0) { $0 + $1.length }
                    + chunk2.reduce(0) { $0 + $1.length }
                if v >= 3 {
                    fputs("[M::main] Processed \(totalPairs) pairs "
                        + "(\(String(format: "%.1f", Double(totalBases) / 1_000_000))M bases)\n",
                        stderr)
                }
                try await aligner.alignPairedBatch(
                    reads1: chunk1, reads2: chunk2,
                    outputFile: outFile, header: header
                )
            }
        } else {
            // Single-end mode: stream from one file
            let inFile = try HTSFile(path: fastqFiles[0], mode: "r")
            let inHeader = try inFile.samHeader()
            let iter = inFile.samIterator(header: inHeader)
            var totalReads = 0
            var totalBases = 0
            if v >= 3 { fputs("Aligning...\n", stderr) }
            while true {
                let chunk = readChunk(
                    from: iter, maxBases: chunkBases,
                    preserveComments: appendComment
                )
                if chunk.isEmpty { break }
                totalReads += chunk.count
                totalBases += chunk.reduce(0) { $0 + $1.length }
                if v >= 3 {
                    fputs("[M::main] Processed \(totalReads) reads "
                        + "(\(String(format: "%.1f", Double(totalBases) / 1_000_000))M bases)\n",
                        stderr)
                }
                try await aligner.alignBatch(
                    reads: chunk, outputFile: outFile, header: header
                )
            }
        }

        if v >= 3 { fputs("Done.\n", stderr) }
    }
}

// MARK: - Streaming FASTQ Reader

/// Default chunk size: ~50M bases per chunk (~330K reads of 150bp, ~100MB memory).
let defaultChunkBases = 50_000_000

/// Parse a BAMRecord into a ReadSequence.
func parseBAMRecord(_ record: borrowing BAMRecord, preserveComments: Bool) -> ReadSequence {
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
    let comment = preserveComments ? (record.auxiliaryData.string(forTag: "CO") ?? "") : ""

    return ReadSequence(name: name, bases: bases, qualities: qualities, comment: comment)
}

/// Read a chunk of reads from an iterator, up to `maxBases` total bases.
func readChunk(
    from iterator: SAMRecordIterator,
    maxBases: Int,
    preserveComments: Bool
) -> [ReadSequence] {
    var reads: [ReadSequence] = []
    var bases = 0
    while let record = iterator.next() {
        let read = parseBAMRecord(record, preserveComments: preserveComments)
        bases += read.length
        reads.append(read)
        if bases >= maxBases { break }
    }
    return reads
}

/// Read exactly `count` records from an iterator (for PE R2 synchronization).
func readExactCount(
    from iterator: SAMRecordIterator,
    count: Int,
    preserveComments: Bool
) -> [ReadSequence] {
    var reads: [ReadSequence] = []
    reads.reserveCapacity(count)
    for _ in 0..<count {
        guard let record = iterator.next() else { break }
        reads.append(parseBAMRecord(record, preserveComments: preserveComments))
    }
    return reads
}

/// Read interleaved pairs (R1, R2, R1, R2, ...) up to `maxBases` bases (counted from R1).
func readInterleavedChunk(
    from iterator: SAMRecordIterator,
    maxBases: Int,
    preserveComments: Bool
) -> ([ReadSequence], [ReadSequence]) {
    var reads1: [ReadSequence] = []
    var reads2: [ReadSequence] = []
    var bases = 0
    while let r1Record = iterator.next() {
        guard let r2Record = iterator.next() else {
            reads1.append(parseBAMRecord(r1Record, preserveComments: preserveComments))
            break
        }
        let r1 = parseBAMRecord(r1Record, preserveComments: preserveComments)
        let r2 = parseBAMRecord(r2Record, preserveComments: preserveComments)
        bases += r1.length
        reads1.append(r1)
        reads2.append(r2)
        if bases >= maxBases { break }
    }
    return (reads1, reads2)
}
