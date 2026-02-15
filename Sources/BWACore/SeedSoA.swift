/// Structure-of-Arrays container for seed data, optimized for hot-path iteration.
///
/// Stores seed fields in separate contiguous columns with 8-byte alignment,
/// enabling cache-efficient sequential access and SIMD vectorization.
///
/// Layout: [rbeg: Int64 × cap][qbeg: Int32 × cap (padded to 8)][len: Int32 × cap (padded to 8)][score: Int32 × cap]
public struct SeedSoA: @unchecked Sendable {
    @usableFromInline let storage: UnsafeMutableRawPointer
    public let capacity: Int
    public private(set) var count: Int

    // Column byte offsets (all 8-byte aligned)
    @usableFromInline let qbegOff: Int
    @usableFromInline let lenOff: Int
    @usableFromInline let scoreOff: Int

    @inlinable public var rbegs: UnsafeMutablePointer<Int64> {
        storage.assumingMemoryBound(to: Int64.self)
    }
    @inlinable public var qbegs: UnsafeMutablePointer<Int32> {
        storage.advanced(by: qbegOff).assumingMemoryBound(to: Int32.self)
    }
    @inlinable public var lens: UnsafeMutablePointer<Int32> {
        storage.advanced(by: lenOff).assumingMemoryBound(to: Int32.self)
    }
    @inlinable public var scores: UnsafeMutablePointer<Int32> {
        storage.advanced(by: scoreOff).assumingMemoryBound(to: Int32.self)
    }

    /// Round up n Int32 elements to 8-byte alignment.
    @usableFromInline
    static func pad4(_ n: Int) -> Int { (n * 4 + 7) & ~7 }

    /// Total bytes needed for a given capacity.
    @inlinable
    static func allocationSize(capacity: Int) -> Int {
        capacity * 8 + pad4(capacity) * 3
    }

    /// Create an empty container with the given capacity.
    public init(capacity: Int) {
        self.capacity = capacity
        self.count = 0
        let size = Self.allocationSize(capacity: capacity)
        self.storage = UnsafeMutableRawPointer.allocate(byteCount: max(size, 1), alignment: 8)
        self.qbegOff = capacity * 8
        self.lenOff = qbegOff + Self.pad4(capacity)
        self.scoreOff = lenOff + Self.pad4(capacity)
    }

    /// Scatter AoS → SoA from an array of MemSeed.
    public init(from seeds: [MemSeed]) {
        self.init(capacity: seeds.count)
        self.count = seeds.count
        let r = rbegs, q = qbegs, l = lens, s = scores
        for i in 0..<seeds.count {
            r[i] = seeds[i].rbeg
            q[i] = seeds[i].qbeg
            l[i] = seeds[i].len
            s[i] = seeds[i].score
        }
    }

    /// Gather SoA → AoS back to an array of MemSeed.
    public func toArray() -> [MemSeed] {
        let r = rbegs, q = qbegs, l = lens, s = scores
        return (0..<count).map { i in
            MemSeed(rbeg: r[i], qbeg: q[i], len: l[i], score: s[i])
        }
    }

    /// Free the backing storage.
    public func deallocate() {
        storage.deallocate()
    }
}
