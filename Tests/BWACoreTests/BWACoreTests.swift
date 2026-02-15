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

    // MARK: - SeedSoA Tests

    @Test("SeedSoA round-trip preserves all fields")
    func testSeedSoARoundTrip() {
        let seeds: [MemSeed] = [
            MemSeed(rbeg: 1000, qbeg: 10, len: 19, score: 19),
            MemSeed(rbeg: 2000, qbeg: 50, len: 25, score: 25),
            MemSeed(rbeg: 3000, qbeg: 80, len: 30, score: 30),
        ]

        let soa = SeedSoA(from: seeds)
        defer { soa.deallocate() }

        #expect(soa.count == 3)
        #expect(soa.capacity == 3)

        let result = soa.toArray()
        #expect(result == seeds)
    }

    @Test("SeedSoA column access matches AoS fields")
    func testSeedSoAColumnAccess() {
        let seeds: [MemSeed] = [
            MemSeed(rbeg: 100, qbeg: 5, len: 20, score: 20),
            MemSeed(rbeg: 200, qbeg: 30, len: 15, score: 15),
        ]

        let soa = SeedSoA(from: seeds)
        defer { soa.deallocate() }

        #expect(soa.rbegs[0] == 100)
        #expect(soa.rbegs[1] == 200)
        #expect(soa.qbegs[0] == 5)
        #expect(soa.qbegs[1] == 30)
        #expect(soa.lens[0] == 20)
        #expect(soa.lens[1] == 15)
        #expect(soa.scores[0] == 20)
        #expect(soa.scores[1] == 15)
    }

    @Test("SeedSoA empty array")
    func testSeedSoAEmpty() {
        let soa = SeedSoA(from: [])
        defer { soa.deallocate() }

        #expect(soa.count == 0)
        #expect(soa.toArray().isEmpty)
    }

    @Test("SeedSoA single element")
    func testSeedSoASingle() {
        let seed = MemSeed(rbeg: 42, qbeg: 7, len: 19, score: 19)
        let soa = SeedSoA(from: [seed])
        defer { soa.deallocate() }

        #expect(soa.count == 1)
        #expect(soa.rbegs[0] == 42)
        #expect(soa.qbegs[0] == 7)
        #expect(soa.lens[0] == 19)
        #expect(soa.scores[0] == 19)
        #expect(soa.toArray() == [seed])
    }

    @Test("SeedSoA large count exercises SIMD path (5+ elements)")
    func testSeedSoALargeCount() {
        let seeds = (0..<13).map { i in
            MemSeed(rbeg: Int64(i) * 100, qbeg: Int32(i) * 10, len: 19, score: 19)
        }

        let soa = SeedSoA(from: seeds)
        defer { soa.deallocate() }

        #expect(soa.count == 13)
        let roundTripped = soa.toArray()
        #expect(roundTripped == seeds)
    }

    @Test("SeedSoA capacity-only init")
    func testSeedSoACapacityInit() {
        let soa = SeedSoA(capacity: 10)
        defer { soa.deallocate() }

        #expect(soa.count == 0)
        #expect(soa.capacity == 10)
    }

    // MARK: - ReadArena Tests

    @Test("ReadArena basic allocation and reset")
    func testReadArenaBasic() {
        var arena = ReadArena(capacity: 4096)
        let ptr1 = arena.allocate(Int64.self, count: 10)
        ptr1[0] = 42
        ptr1[9] = 99
        #expect(ptr1[0] == 42)
        #expect(ptr1[9] == 99)

        let ptr2 = arena.allocate(UInt8.self, count: 100)
        ptr2[0] = 7
        #expect(ptr2[0] == 7)

        // Reset reuses the same buffer
        arena.reset()
        let ptr3 = arena.allocate(Int64.self, count: 10)
        ptr3[0] = 123
        #expect(ptr3[0] == 123)
    }

    @Test("ReadArena alignment")
    func testReadArenaAlignment() {
        var arena = ReadArena(capacity: 4096)
        // Allocate a byte to offset, then allocate Int64 which needs 8-byte alignment
        let _ = arena.allocate(UInt8.self, count: 1)
        let ptr = arena.allocate(Int64.self, count: 1)
        #expect(Int(bitPattern: ptr) % MemoryLayout<Int64>.alignment == 0)
    }

    // MARK: - ArenaBuffer Tests

    @Test("ArenaBuffer append and subscript")
    func testArenaBufferAppendSubscript() {
        var arena = ReadArena(capacity: 4096)
        var buf = ArenaBuffer<Int32>(
            base: arena.allocate(Int32.self, count: 10),
            capacity: 10
        )
        #expect(buf.count == 0)

        buf.append(10)
        buf.append(20)
        buf.append(30)
        #expect(buf.count == 3)
        #expect(buf[0] == 10)
        #expect(buf[1] == 20)
        #expect(buf[2] == 30)

        buf[1] = 99
        #expect(buf[1] == 99)
    }

    @Test("ArenaBuffer removeAll resets count")
    func testArenaBufferRemoveAll() {
        var arena = ReadArena(capacity: 4096)
        var buf = ArenaBuffer<Int32>(
            base: arena.allocate(Int32.self, count: 10),
            capacity: 10
        )
        buf.append(1)
        buf.append(2)
        #expect(buf.count == 2)

        buf.removeAll()
        #expect(buf.count == 0)

        // Can reuse capacity
        buf.append(3)
        #expect(buf.count == 1)
        #expect(buf[0] == 3)
    }

    @Test("ArenaBuffer sort")
    func testArenaBufferSort() {
        var arena = ReadArena(capacity: 4096)
        var buf = ArenaBuffer<Int32>(
            base: arena.allocate(Int32.self, count: 10),
            capacity: 10
        )
        buf.append(30)
        buf.append(10)
        buf.append(20)

        buf.sort(by: <)
        #expect(buf[0] == 10)
        #expect(buf[1] == 20)
        #expect(buf[2] == 30)
    }

    @Test("ArenaBuffer with SMEM type")
    func testArenaBufferSMEM() {
        var arena = ReadArena(capacity: 4096)
        var buf = ArenaBuffer<SMEM>(
            base: arena.allocate(SMEM.self, count: 10),
            capacity: 10
        )

        buf.append(SMEM(k: 100, l: 200, queryBegin: 5, queryEnd: 15))
        buf.append(SMEM(k: 300, l: 400, queryBegin: 0, queryEnd: 10))
        #expect(buf.count == 2)

        // Sort by queryBegin
        buf.sort { $0.queryBegin < $1.queryBegin }
        #expect(buf[0].queryBegin == 0)
        #expect(buf[1].queryBegin == 5)

        // Convert to Array
        let arr = Array(UnsafeBufferPointer(start: buf.storage, count: buf.count))
        #expect(arr.count == 2)
        #expect(arr[0].k == 300)
    }
}
