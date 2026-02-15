/// 2-bit nucleotide encoding matching BWA convention: A=0, C=1, G=2, T=3, N=4
public enum Base: UInt8, Sendable {
    case A = 0
    case C = 1
    case G = 2
    case T = 3
    case N = 4

    /// Complement base (A<->T, C<->G)
    @inlinable
    public var complement: Base {
        switch self {
        case .A: return .T
        case .T: return .A
        case .C: return .G
        case .G: return .C
        case .N: return .N
        }
    }

    /// Decode from ASCII character
    @inlinable
    public init(ascii: UInt8) {
        switch ascii {
        case UInt8(ascii: "A"), UInt8(ascii: "a"): self = .A
        case UInt8(ascii: "C"), UInt8(ascii: "c"): self = .C
        case UInt8(ascii: "G"), UInt8(ascii: "g"): self = .G
        case UInt8(ascii: "T"), UInt8(ascii: "t"): self = .T
        default: self = .N
        }
    }

}

/// Lookup table for ASCII -> 2-bit encoding (matches bwa's nst_nt4_table)
/// Returns 4 for non-ACGT characters
public let nst_nt4_table: [UInt8] = {
    var table = [UInt8](repeating: 4, count: 256)
    table[Int(UInt8(ascii: "A"))] = 0; table[Int(UInt8(ascii: "a"))] = 0
    table[Int(UInt8(ascii: "C"))] = 1; table[Int(UInt8(ascii: "c"))] = 1
    table[Int(UInt8(ascii: "G"))] = 2; table[Int(UInt8(ascii: "g"))] = 2
    table[Int(UInt8(ascii: "T"))] = 3; table[Int(UInt8(ascii: "t"))] = 3
    return table
}()
