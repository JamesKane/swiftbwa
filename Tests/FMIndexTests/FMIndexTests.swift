import Testing
@testable import FMIndex
@testable import BWACore

@Suite("FMIndex Tests")
struct FMIndexTests {

    @Test("CheckpointOCC size is 64 bytes")
    func testCheckpointSize() {
        // CP_OCC should be 64 bytes: 4x Int64 (32) + 4x UInt64 (32)
        #expect(MemoryLayout<CheckpointOCC>.size == 64)
    }

    @Test("CP constants are correct")
    func testConstants() {
        #expect(CP_BLOCK_SIZE == 64)
        #expect(CP_SHIFT == 6)
        #expect(CP_MASK == 63)
    }

    @Test("BWT one-hot mask array construction")
    func testOneHotMaskArray() {
        // Create a minimal BWT to test mask construction
        let cpBuffer = UnsafeMutablePointer<CheckpointOCC>.allocate(capacity: 1)
        cpBuffer.initialize(to: CheckpointOCC())
        let bwt = BWT(
            checkpoints: UnsafeBufferPointer(start: cpBuffer, count: 1),
            count: (1, 2, 3, 4, 5),
            length: 64,
            sentinelIndex: 0,
            ownedBase: UnsafeMutableRawPointer(cpBuffer)
        )

        #expect(bwt.oneHotMaskArray[0] == 0)
        #expect(bwt.oneHotMaskArray[1] == 0x8000_0000_0000_0000)
        #expect(bwt.oneHotMaskArray[2] == 0xC000_0000_0000_0000)
        // mask[63] has bits 1..63 set (63 bits), bit 0 is clear
        #expect(bwt.oneHotMaskArray[63] == UInt64.max - 1)
    }

    @Test("BWT count tuple accessors")
    func testBWTCountAccessors() {
        let cpBuffer = UnsafeMutablePointer<CheckpointOCC>.allocate(capacity: 1)
        cpBuffer.initialize(to: CheckpointOCC())
        let bwt = BWT(
            checkpoints: UnsafeBufferPointer(start: cpBuffer, count: 1),
            count: (10, 20, 30, 40, 50),
            length: 64,
            sentinelIndex: 0,
            ownedBase: UnsafeMutableRawPointer(cpBuffer)
        )

        #expect(bwt.count(for: 0) == 10)
        #expect(bwt.count(for: 1) == 20)
        #expect(bwt.count(for: 2) == 30)
        #expect(bwt.count(for: 3) == 40)
        #expect(bwt.count(for: 4) == 50)
        #expect(bwt.count(forNext: 0) == 20)
        #expect(bwt.count(forNext: 3) == 50)
    }

    @Test("SuffixArray entry reconstruction")
    func testSuffixArrayEntry() {
        let msBuffer = UnsafeMutablePointer<Int8>.allocate(capacity: 2)
        let lsBuffer = UnsafeMutablePointer<UInt32>.allocate(capacity: 2)
        msBuffer[0] = 0; lsBuffer[0] = 100      // SA[0] = 100
        msBuffer[1] = 1; lsBuffer[1] = 500      // SA[1] = (1 << 32) + 500

        let sa = SuffixArray(
            msBytes: UnsafeBufferPointer(start: msBuffer, count: 2),
            lsWords: UnsafeBufferPointer(start: lsBuffer, count: 2),
            count: 2,
            compressionShift: 0,
            ownedMSBase: UnsafeMutableRawPointer(msBuffer),
            ownedLSBase: UnsafeMutableRawPointer(lsBuffer)
        )

        #expect(sa.entry(at: 0) == 100)
        #expect(sa.entry(at: 1) == (1 << 32) + 500)
    }

    @Test("PackedReference base extraction")
    func testPackedReference() {
        // Pack "ACGT" = 0b00_01_10_11 = 0x1B
        let buffer = UnsafeMutablePointer<UInt8>.allocate(capacity: 1)
        buffer[0] = 0b00_01_10_11  // A=00, C=01, G=10, T=11

        let pac = PackedReference(
            data: UnsafeBufferPointer(start: buffer, count: 1),
            length: 4,
            ownedBase: UnsafeMutableRawPointer(buffer)
        )

        #expect(pac.base(at: 0) == 0)  // A
        #expect(pac.base(at: 1) == 1)  // C
        #expect(pac.base(at: 2) == 2)  // G
        #expect(pac.base(at: 3) == 3)  // T
    }

    @Test("PackedReference subsequence")
    func testPackedSubsequence() {
        // Pack "ACGTACGT" = 2 bytes: 0x1B, 0x1B
        let buffer = UnsafeMutablePointer<UInt8>.allocate(capacity: 2)
        buffer[0] = 0b00_01_10_11
        buffer[1] = 0b00_01_10_11

        let pac = PackedReference(
            data: UnsafeBufferPointer(start: buffer, count: 2),
            length: 8,
            ownedBase: UnsafeMutableRawPointer(buffer)
        )

        let sub = pac.subsequence(from: 2, length: 4)
        #expect(sub == [2, 3, 0, 1])  // GTAC
    }

    @Test("BackwardSearch initInterval")
    func testInitInterval() {
        // Create a BWT with known counts: A=10..20, C=20..30, G=30..40, T=40..50
        let cpBuffer = UnsafeMutablePointer<CheckpointOCC>.allocate(capacity: 1)
        cpBuffer.initialize(to: CheckpointOCC())
        let bwt = BWT(
            checkpoints: UnsafeBufferPointer(start: cpBuffer, count: 1),
            count: (10, 20, 30, 40, 50),
            length: 64,
            sentinelIndex: 0,
            ownedBase: UnsafeMutableRawPointer(cpBuffer)
        )

        let intA = BackwardSearch.initInterval(bwt: bwt, base: 0)
        #expect(intA.k == 10)
        #expect(intA.s == 10)

        let intC = BackwardSearch.initInterval(bwt: bwt, base: 1)
        #expect(intC.k == 20)
        #expect(intC.s == 10)

        let intN = BackwardSearch.initInterval(bwt: bwt, base: 5)
        #expect(intN.s == 0)
    }
}
