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
}
