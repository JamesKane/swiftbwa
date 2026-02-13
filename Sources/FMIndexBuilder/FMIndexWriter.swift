import Foundation
import BWACore
import FMIndex

/// Derives BWT from suffix array, builds FM-index checkpoints and compressed SA,
/// and writes the .bwt.2bit.64 file matching bwa-mem2's format exactly.
struct FMIndexWriter: Sendable {

    private static let SA_COMPX = 3  // sample every 8th SA entry
    private static let DUMMY_CHAR: UInt8 = 6  // padding character for 64-aligned BWT

    /// Write the .bwt.2bit.64 file from the suffix array and binary sequence.
    ///
    /// - Parameters:
    ///   - prefix: Output file prefix (writes `prefix.bwt.2bit.64`)
    ///   - binarySeq: The full forward+reverse-complement sequence (0-3 encoded, length = 2*pacLen)
    ///   - pacLen: Forward strand length
    ///   - suffixArray: SA of length (2*pacLen + 2), where sa[0] = 2*pacLen (sentinel)
    ///   - counts: Character frequency counts [A, C, G, T, $] for the full BWT string
    ///   - verbosity: Logging level
    static func write(
        prefix: String,
        binarySeq: UnsafePointer<UInt8>,
        pacLen: Int64,
        suffixArray: UnsafeMutablePointer<Int64>,
        counts: inout [Int64],
        verbosity: Int
    ) throws {
        let refSeqLen = 2 * pacLen + 1  // fwd + rc + sentinel
        let n = refSeqLen

        // Step 1: Derive BWT from SA
        // bwt[i] = binarySeq[sa[i]-1] if sa[i]>0, else sentinel(=4)
        // Pad to 64-aligned length with DUMMY_CHAR
        let bwtAlignedLen = ((n + 63) / 64) * 64
        let bwt = UnsafeMutablePointer<UInt8>.allocate(capacity: Int(bwtAlignedLen))
        defer { bwt.deallocate() }

        var sentinelIndex: Int64 = -1
        for i in 0..<Int(n) {
            let saVal = suffixArray[i]
            if saVal == 0 {
                bwt[i] = 4  // sentinel
                sentinelIndex = Int64(i)
            } else {
                bwt[i] = binarySeq[Int(saVal - 1)]
            }
        }
        // Pad remainder with DUMMY_CHAR
        for i in Int(n)..<Int(bwtAlignedLen) {
            bwt[i] = DUMMY_CHAR
        }

        if verbosity >= 3 {
            fputs("[M::index] BWT constructed, sentinel at position \(sentinelIndex)\n", stderr)
        }

        // Step 2: Build cumulative counts (before +1 adjustment)
        // counts[0..4] = [A, C, G, T, $] frequencies
        // For the file, we write cumulative counts: count[i] = sum of counts[0..<i]
        var cumCounts = [Int64](repeating: 0, count: 5)
        var sum: Int64 = 0
        for i in 0..<5 {
            cumCounts[i] = sum
            sum += counts[i]
        }

        // Step 3: Build checkpoints
        let cpCount = Int(n >> CP_SHIFT) + 1
        let checkpoints = UnsafeMutablePointer<CheckpointOCC>.allocate(capacity: cpCount)
        defer { checkpoints.deallocate() }

        var runCounts: [Int64] = [0, 0, 0, 0]  // running A, C, G, T counts
        for cpIdx in 0..<cpCount {
            var cp = CheckpointOCC()
            cp.counts = (runCounts[0], runCounts[1], runCounts[2], runCounts[3])

            // Build one-hot bitstrings for this 64-char block
            var bs0: UInt64 = 0  // A
            var bs1: UInt64 = 0  // C
            var bs2: UInt64 = 0  // G
            var bs3: UInt64 = 0  // T

            let blockStart = cpIdx * CP_BLOCK_SIZE
            for pos in 0..<CP_BLOCK_SIZE {
                let globalPos = blockStart + pos
                guard globalPos < Int(n) else { break }
                let ch = bwt[globalPos]
                // MSB = first base in block: bit (63 - pos) set
                let bit: UInt64 = 1 << (63 - pos)
                switch ch {
                case 0: bs0 |= bit; runCounts[0] += 1
                case 1: bs1 |= bit; runCounts[1] += 1
                case 2: bs2 |= bit; runCounts[2] += 1
                case 3: bs3 |= bit; runCounts[3] += 1
                default: break  // sentinel or padding
                }
            }

            cp.bitstrings = (bs0, bs1, bs2, bs3)
            checkpoints[cpIdx] = cp
        }

        if verbosity >= 3 {
            fputs("[M::index] Built \(cpCount) checkpoints\n", stderr)
        }

        // Step 4: Compress SA (SA_COMPX=3: sample every 8th entry)
        let saCount = Int(n >> SA_COMPX) + 1
        let msByte = UnsafeMutablePointer<Int8>.allocate(capacity: saCount)
        let lsWord = UnsafeMutablePointer<UInt32>.allocate(capacity: saCount)
        defer { msByte.deallocate(); lsWord.deallocate() }

        for i in 0..<saCount {
            let saIdx = i << SA_COMPX
            let saVal: Int64
            if saIdx < Int(n) {
                saVal = suffixArray[saIdx]
            } else {
                saVal = 0
            }
            msByte[i] = Int8(truncatingIfNeeded: (saVal >> 32) & 0xFF)
            lsWord[i] = UInt32(truncatingIfNeeded: saVal & 0xFFFF_FFFF)
        }

        if verbosity >= 3 {
            fputs("[M::index] Compressed SA: \(saCount) entries\n", stderr)
        }

        // Step 5: Write .bwt.2bit.64 file
        let bwtPath = prefix + ".bwt.2bit.64"
        guard let file = fopen(bwtPath, "wb") else {
            throw BWAError.fileIOError("Cannot create: \(bwtPath)")
        }
        defer { fclose(file) }

        // Header: ref_seq_len (Int64)
        var refLen = refSeqLen
        fwrite(&refLen, 8, 1, file)

        // Cumulative counts: 5 × Int64
        for i in 0..<5 {
            var val = cumCounts[i]
            fwrite(&val, 8, 1, file)
        }

        // Checkpoints
        fwrite(checkpoints, MemoryLayout<CheckpointOCC>.size, cpCount, file)

        // SA ms_byte
        fwrite(msByte, 1, saCount, file)

        // SA ls_word
        fwrite(lsWord, 4, saCount, file)

        // Sentinel index (Int64)
        var sentinel = sentinelIndex
        fwrite(&sentinel, 8, 1, file)

        if verbosity >= 3 {
            fputs("[M::index] Wrote \(bwtPath)\n", stderr)
        }
    }
}
