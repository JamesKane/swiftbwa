/// A sequencing read with 2-bit encoded bases and quality scores.
public struct ReadSequence: Sendable {
    /// Read name
    public var name: String
    /// 2-bit encoded sequence (A=0, C=1, G=2, T=3, N=4)
    public var bases: [UInt8]
    /// Phred quality scores
    public var qualities: [UInt8]
    /// Optional comment/description from FASTQ header
    public var comment: String

    /// Length of the read
    @inlinable
    public var length: Int { bases.count }

    public init(name: String, bases: [UInt8], qualities: [UInt8], comment: String = "") {
        self.name = name
        self.bases = bases
        self.qualities = qualities
        self.comment = comment
    }

    /// Create from ASCII FASTQ sequence and quality strings
    public init(name: String, sequence: String, qualityString: String, comment: String = "") {
        self.name = name
        self.comment = comment
        self.bases = sequence.utf8.map { nst_nt4_table[Int($0)] }
        self.qualities = qualityString.utf8.map { $0 - 33 }  // Phred+33 -> raw
    }

    /// Reverse complement of this read
    public func reverseComplement() -> [UInt8] {
        bases.reversed().map { b in b < 4 ? 3 - b : 4 }
    }
}
