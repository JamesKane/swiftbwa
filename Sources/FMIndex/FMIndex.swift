import BWACore

/// Combined FM-Index containing BWT, suffix array, packed reference, and metadata.
public final class FMIndex: @unchecked Sendable {
    public let bwt: BWT
    public let suffixArray: SuffixArray
    public let packedRef: PackedReference
    public let metadata: ReferenceMetadata

    public init(bwt: BWT, suffixArray: SuffixArray, packedRef: PackedReference, metadata: ReferenceMetadata) {
        self.bwt = bwt
        self.suffixArray = suffixArray
        self.packedRef = packedRef
        self.metadata = metadata
    }

    /// Reference sequence length (2 * genome_length + 1, the BWT length)
    public var referenceSeqLen: Int64 { bwt.length }

    /// Genome length (one strand only, = referenceSeqLen / 2)
    public var genomeLength: Int64 { bwt.length / 2 }

    /// Fetch reference bases into a pre-allocated buffer, handling both strands.
    /// Returns the number of bases written.
    public func getReference(at pos: Int64, length: Int, into buffer: UnsafeMutablePointer<UInt8>) -> Int {
        let gl = genomeLength
        if pos >= gl {
            let fwdEnd = 2 * gl - pos
            let fwdStart = fwdEnd - Int64(length)
            let safeStart = max(0, fwdStart)
            let safeEnd = min(fwdEnd, gl)
            let safeLen = Int(safeEnd - safeStart)
            guard safeLen > 0 else { return 0 }
            packedRef.subsequence(from: safeStart, length: safeLen, into: buffer)
            var lo = 0; var hi = safeLen - 1
            while lo < hi {
                let tmp = buffer[lo]; buffer[lo] = buffer[hi]; buffer[hi] = tmp
                lo += 1; hi -= 1
            }
            for i in 0..<safeLen { buffer[i] = 3 - buffer[i] }
            return safeLen
        } else {
            let safeLen = min(length, Int(gl - pos))
            guard safeLen > 0 && pos >= 0 else { return 0 }
            packedRef.subsequence(from: pos, length: safeLen, into: buffer)
            return safeLen
        }
    }

    /// Fetch reference bases at a BWT-space position, handling both strands.
    ///
    /// The packed reference (.pac) stores only the forward strand. Positions in
    /// `[0, genomeLength)` access the forward strand directly. Positions in
    /// `[genomeLength, 2*genomeLength)` are on the reverse complement strand:
    /// the corresponding forward bases are extracted, reversed, and complemented.
    ///
    /// Matches bwa-mem2's `bns_fetch_seq()` / `bns_get_seq()` behavior.
    public func getReference(at pos: Int64, length: Int) -> [UInt8] {
        let gl = genomeLength
        if pos >= gl {
            // Reverse strand: BWT rev range [pos, pos+length)
            // maps to forward range [2*gl - pos - length, 2*gl - pos)
            let fwdEnd = 2 * gl - pos
            let fwdStart = fwdEnd - Int64(length)
            let safeStart = max(0, fwdStart)
            let safeEnd = min(fwdEnd, gl)
            let safeLen = Int(safeEnd - safeStart)
            guard safeLen > 0 else { return [] }
            var bases = packedRef.subsequence(from: safeStart, length: safeLen)
            bases.reverse()
            for i in 0..<bases.count {
                bases[i] = 3 - bases[i]  // A↔T, C↔G
            }
            return bases
        } else {
            let safeLen = min(length, Int(gl - pos))
            guard safeLen > 0 && pos >= 0 else { return [] }
            return packedRef.subsequence(from: pos, length: safeLen)
        }
    }
}
