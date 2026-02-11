import Testing
@testable import BWACore

@Suite("BWACore Tests")
struct BWACoreTests {

    @Test("Base encoding from ASCII")
    func testBaseEncoding() {
        #expect(Base(ascii: UInt8(ascii: "A")) == .A)
        #expect(Base(ascii: UInt8(ascii: "C")) == .C)
        #expect(Base(ascii: UInt8(ascii: "G")) == .G)
        #expect(Base(ascii: UInt8(ascii: "T")) == .T)
        #expect(Base(ascii: UInt8(ascii: "N")) == .N)
        #expect(Base(ascii: UInt8(ascii: "a")) == .A)
        #expect(Base(ascii: UInt8(ascii: "X")) == .N)
    }

    @Test("Base complement")
    func testBaseComplement() {
        #expect(Base.A.complement == .T)
        #expect(Base.T.complement == .A)
        #expect(Base.C.complement == .G)
        #expect(Base.G.complement == .C)
        #expect(Base.N.complement == .N)
    }

    @Test("nst_nt4_table encoding")
    func testNstNt4Table() {
        #expect(nst_nt4_table[Int(UInt8(ascii: "A"))] == 0)
        #expect(nst_nt4_table[Int(UInt8(ascii: "C"))] == 1)
        #expect(nst_nt4_table[Int(UInt8(ascii: "G"))] == 2)
        #expect(nst_nt4_table[Int(UInt8(ascii: "T"))] == 3)
        #expect(nst_nt4_table[Int(UInt8(ascii: "N"))] == 4)
        #expect(nst_nt4_table[Int(UInt8(ascii: "Z"))] == 4)
    }

    @Test("ScoringParameters defaults match BWA-MEM")
    func testScoringDefaults() {
        let params = ScoringParameters()
        #expect(params.matchScore == 1)
        #expect(params.mismatchPenalty == 4)
        #expect(params.gapOpenPenalty == 6)
        #expect(params.gapExtendPenalty == 1)
        #expect(params.bandWidth == 100)
        #expect(params.zDrop == 100)
        #expect(params.minSeedLength == 19)
        #expect(params.maxOccurrences == 500)
        #expect(params.minOutputScore == 30)
        #expect(params.maxChainGap == 10000)
    }

    @Test("Scoring matrix construction")
    func testScoringMatrix() {
        let params = ScoringParameters()
        let mat = params.scoringMatrix()
        #expect(mat.count == 25)
        // Diagonal: match score = 1
        #expect(mat[0] == 1)   // A-A
        #expect(mat[6] == 1)   // C-C
        #expect(mat[12] == 1)  // G-G
        #expect(mat[18] == 1)  // T-T
        // Off-diagonal: -mismatch = -4
        #expect(mat[1] == -4)  // A-C
        #expect(mat[5] == -4)  // C-A
        // N column/row: -1
        #expect(mat[4] == -1)  // A-N
        #expect(mat[20] == -1) // N-A
        #expect(mat[24] == -1) // N-N
    }

    @Test("ReadSequence from ASCII strings")
    func testReadSequence() {
        let read = ReadSequence(name: "test", sequence: "ACGT", qualityString: "IIII")
        #expect(read.length == 4)
        #expect(read.bases == [0, 1, 2, 3])
        #expect(read.qualities == [40, 40, 40, 40])  // 'I' - 33 = 40
    }

    @Test("ReadSequence reverse complement")
    func testReverseComplement() {
        let read = ReadSequence(name: "test", sequence: "ACGT", qualityString: "IIII")
        let rc = read.reverseComplement()
        #expect(rc == [0, 1, 2, 3])  // ACGT -> ACGT (palindrome)

        let read2 = ReadSequence(name: "test", sequence: "AACG", qualityString: "IIII")
        let rc2 = read2.reverseComplement()
        #expect(rc2 == [1, 2, 3, 3])  // AACG -> CGTT
    }

    @Test("SMEM properties")
    func testSMEM() {
        let smem = SMEM(k: 10, l: 20, queryBegin: 5, queryEnd: 15)
        #expect(smem.count == 11)
        #expect(smem.length == 10)
    }

    @Test("ReferenceMetadata sequenceID binary search")
    func testSequenceID() {
        var metadata = ReferenceMetadata()
        metadata.numSequences = 3
        metadata.annotations = [
            ReferenceAnnotation(offset: 0, length: 100, name: "chr1"),
            ReferenceAnnotation(offset: 100, length: 200, name: "chr2"),
            ReferenceAnnotation(offset: 300, length: 150, name: "chr3"),
        ]

        #expect(metadata.sequenceID(for: 0) == 0)
        #expect(metadata.sequenceID(for: 99) == 0)
        #expect(metadata.sequenceID(for: 100) == 1)
        #expect(metadata.sequenceID(for: 299) == 1)
        #expect(metadata.sequenceID(for: 300) == 2)
        #expect(metadata.sequenceID(for: 449) == 2)
    }

    @Test("ReferenceMetadata decodePosition")
    func testDecodePosition() {
        var metadata = ReferenceMetadata()
        metadata.numSequences = 2
        metadata.annotations = [
            ReferenceAnnotation(offset: 0, length: 100, name: "chr1"),
            ReferenceAnnotation(offset: 100, length: 200, name: "chr2"),
        ]

        let (rid, pos) = metadata.decodePosition(150)
        #expect(rid == 1)
        #expect(pos == 50)
    }
}
