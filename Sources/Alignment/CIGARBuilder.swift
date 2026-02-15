import BWACore

/// CIGAR operation codes matching BAM spec.
public enum CIGAROp: UInt32, Sendable {
    case match = 0      // M
    case insertion = 1   // I
    case deletion = 2    // D
    case softClip = 4    // S
    case hardClip = 5    // H
}

/// Builds a CIGAR array from alignment operations.
public struct CIGARBuilder: Sendable {
    private var operations: [(op: CIGAROp, length: UInt32)]

    public init() {
        operations = []
    }

    /// Append an operation. Merges with previous if same op.
    public mutating func append(_ op: CIGAROp, length: UInt32 = 1) {
        if let last = operations.last, last.op == op {
            operations[operations.count - 1].length += length
        } else {
            operations.append((op, length))
        }
    }

    /// Build packed CIGAR array (each UInt32 = length << 4 | op)
    public func build() -> [UInt32] {
        operations.map { UInt32($0.length) << 4 | $0.op.rawValue }
    }

    /// Reverse the order of operations (for left extension)
    public mutating func reverse() {
        operations.reverse()
    }
}
