/// Per-read arena allocator. One 64KB malloc replaces ~150+ per-read mallocs.
/// `~Copyable` ensures unique ownership and deterministic deallocation.
public struct ReadArena: ~Copyable {
    @usableFromInline let buffer: UnsafeMutableRawPointer
    @usableFromInline let totalSize: Int
    @usableFromInline var offset: Int = 0

    @inlinable
    public init(capacity: Int = 65536) {
        self.totalSize = capacity
        self.buffer = .allocate(byteCount: capacity, alignment: 16)
    }

    /// Bump-allocate a typed buffer. ~2 instructions vs malloc's ~100+.
    @inline(__always)
    public mutating func allocate<T>(_ type: T.Type, count: Int) -> UnsafeMutablePointer<T> {
        let alignment = MemoryLayout<T>.alignment
        let aligned = (offset + alignment - 1) & ~(alignment - 1)
        let bytes = MemoryLayout<T>.stride * count
        precondition(aligned + bytes <= totalSize, "ReadArena overflow")
        let ptr = buffer.advanced(by: aligned).bindMemory(to: T.self, capacity: count)
        offset = aligned + bytes
        return ptr
    }

    /// Reset for next read — zero cost, just resets the bump pointer.
    @inline(__always)
    public mutating func reset() { offset = 0 }

    deinit { buffer.deallocate() }
}

/// Non-owning typed buffer backed by arena memory. Supports append, subscript, sort.
public struct ArenaBuffer<Element> {
    public let storage: UnsafeMutablePointer<Element>
    public let capacity: Int
    public private(set) var count: Int = 0

    @inlinable
    public init(base: UnsafeMutablePointer<Element>, capacity: Int) {
        self.storage = base
        self.capacity = capacity
    }

    @inline(__always)
    public mutating func append(_ element: Element) {
        precondition(count < capacity, "ArenaBuffer overflow")
        storage[count] = element
        count &+= 1
    }

    @inline(__always)
    public mutating func removeAll() { count = 0 }

    public subscript(index: Int) -> Element {
        @inline(__always) get { storage[index] }
        @inline(__always) set { storage[index] = newValue }
    }

    /// Sort the filled portion in-place.
    @inlinable
    public mutating func sort(by areInIncreasingOrder: (Element, Element) -> Bool) {
        var buf = UnsafeMutableBufferPointer(start: storage, count: count)
        buf.sort(by: areInIncreasingOrder)
    }
}
