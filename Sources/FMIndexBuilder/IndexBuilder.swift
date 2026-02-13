import Foundation
import BWACore

/// Top-level orchestrator for building a BWA-MEM2 compatible index from a FASTA file.
///
/// Produces 4 files: `.pac`, `.ann`, `.amb`, `.bwt.2bit.64`
public struct IndexBuilder: Sendable {

    /// Build a complete BWA-MEM2 index from a FASTA file.
    ///
    /// - Parameters:
    ///   - fasta: Path to the input FASTA file (plain or gzipped)
    ///   - prefix: Output file prefix (files will be `prefix.pac`, etc.)
    ///   - verbosity: Logging verbosity (3=info)
    public static func build(fasta: String, prefix: String, verbosity: Int) throws {
        let startTime = DispatchTime.now()

        // Stage 1: FASTA → .pac + .ann + .amb
        if verbosity >= 3 {
            fputs("[M::index] Stage 1: Encoding FASTA to .pac/.ann/.amb\n", stderr)
        }
        let result = try FASTAEncoder.encode(fasta: fasta, prefix: prefix, verbosity: verbosity)
        let pacLen = result.pacLen
        if verbosity >= 3 {
            fputs("[M::index] Forward strand: \(pacLen) bases, "
                  + "\(result.metadata.numSequences) sequences, "
                  + "\(result.metadata.ambiguities.count) ambiguity regions\n", stderr)
        }

        // Stage 2: Decode .pac → binary_seq (fwd) + append reverse complement
        if verbosity >= 3 {
            fputs("[M::index] Stage 2: Building forward+reverse-complement sequence\n", stderr)
        }
        let seqLen = 2 * pacLen  // fwd + rc (sentinel added by SA construction)
        let binarySeq = UnsafeMutablePointer<UInt8>.allocate(capacity: Int(seqLen))
        defer { binarySeq.deallocate() }

        // Read pac file back to decode
        let pacPath = prefix + ".pac"
        guard let pacFile = fopen(pacPath, "rb") else {
            throw BWAError.fileIOError("Cannot open: \(pacPath)")
        }
        defer { fclose(pacFile) }

        // Read pac bytes (file size - 1 for remainder byte, possibly - 1 more for zero pad)
        fseek(pacFile, 0, SEEK_END)
        let pacFileSize = ftell(pacFile)
        fseek(pacFile, 0, SEEK_SET)

        let pacBytes = UnsafeMutablePointer<UInt8>.allocate(capacity: pacFileSize)
        defer { pacBytes.deallocate() }
        fread(pacBytes, 1, pacFileSize, pacFile)

        // Decode 2-bit packed bases to byte array
        for i in 0..<Int(pacLen) {
            let byteIdx = i >> 2
            let shift = ((~i) & 3) << 1
            binarySeq[i] = (pacBytes[byteIdx] >> shift) & 3
        }

        // Append reverse complement
        // RC of base b (0=A,1=C,2=G,3=T) is 3-b
        for i in 0..<Int(pacLen) {
            binarySeq[Int(pacLen) + i] = 3 - binarySeq[Int(pacLen) - 1 - i]
        }

        // Stage 3: Build suffix array via SAIS
        if verbosity >= 3 {
            fputs("[M::index] Stage 3: Building suffix array (\(seqLen + 1) entries)\n", stderr)
        }

        // The text for SA construction is binarySeq[0..<seqLen] with a sentinel appended.
        // bwa-mem2 convention: SA has seqLen+1 entries, sa[0] is reserved for the sentinel position.
        // We construct SA for text of length seqLen+1 where text[seqLen] = sentinel (encoded as a value
        // smaller than any base, which we handle by setting it to 0 and using alphabetSize that includes it).
        //
        // Actually, bwa-mem2 uses: sa[0]=seqLen, then SAIS on binarySeq[0..<seqLen] with sa[1..].
        // The final SA has seqLen+1 entries.
        let saLen = Int(seqLen) + 1
        let suffixArray = UnsafeMutablePointer<Int64>.allocate(capacity: saLen + 1)
        defer { suffixArray.deallocate() }

        // sa[0] = seqLen (sentinel position)
        suffixArray[0] = seqLen

        // Count character frequencies for the BWT
        // counts[0..3] = A, C, G, T frequencies in binarySeq
        var counts = [Int64](repeating: 0, count: 5)
        for i in 0..<Int(seqLen) {
            counts[Int(binarySeq[i])] += 1
        }
        counts[4] = 1  // sentinel

        // Run SAIS on binarySeq, storing result in sa[1..]
        let saSlice = suffixArray + 1
        let saResult = UnsafeMutableBufferPointer(start: saSlice, count: Int(seqLen))
        binarySeq.withMemoryRebound(to: UInt8.self, capacity: Int(seqLen)) { textPtr in
            saResult.baseAddress!.withMemoryRebound(to: Int64.self, capacity: Int(seqLen)) { saPtr in
                // Zero initialize
                for i in 0..<Int(seqLen) { saPtr[i] = 0 }
            }
            // SAIS builds SA for text[0..<seqLen] with alphabet {0,1,2,3}
            let sa = SAIS.buildSuffixArray(text: textPtr, n: seqLen, alphabetSize: 4)
            for i in 0..<Int(seqLen) {
                saSlice[i] = sa[i]
            }
        }

        if verbosity >= 3 {
            fputs("[M::index] Suffix array constructed\n", stderr)
        }

        // Stage 4: BWT + checkpoints + compressed SA → .bwt.2bit.64
        if verbosity >= 3 {
            fputs("[M::index] Stage 4: Writing .bwt.2bit.64\n", stderr)
        }
        try FMIndexWriter.write(
            prefix: prefix,
            binarySeq: binarySeq,
            pacLen: pacLen,
            suffixArray: suffixArray,
            counts: &counts,
            verbosity: verbosity
        )

        let elapsed = Double(DispatchTime.now().uptimeNanoseconds - startTime.uptimeNanoseconds) / 1e9
        if verbosity >= 3 {
            fputs("[M::index] Index built in \(String(format: "%.1f", elapsed)) seconds\n", stderr)
        }
    }
}
