/// Annotation for a single reference sequence (from .ann file).
///
/// Corresponds to `bntann1_t` in bwa's `bntseq.h`.
public struct ReferenceAnnotation: Sendable {
    /// Offset in the concatenated reference
    public var offset: Int64
    /// Length of this sequence
    public var length: Int32
    /// Number of ambiguous bases (Ns)
    public var nAmb: Int32
    /// GI number
    public var gi: UInt32
    /// Sequence name
    public var name: String
    /// Annotation string (e.g., description line)
    public var anno: String

    public init(offset: Int64 = 0, length: Int32 = 0, nAmb: Int32 = 0,
                gi: UInt32 = 0, name: String = "", anno: String = "") {
        self.offset = offset
        self.length = length
        self.nAmb = nAmb
        self.gi = gi
        self.name = name
        self.anno = anno
    }
}

/// Ambiguity region (stretch of Ns) in the reference (from .amb file).
///
/// Corresponds to `bntamb1_t` in bwa's `bntseq.h`.
public struct AmbiguityRegion: Sendable {
    /// Offset in the concatenated reference
    public var offset: Int64
    /// Length of the ambiguous region
    public var length: Int32
    /// The ambiguous character (encoded as 2-bit)
    public var amb: UInt8

    public init(offset: Int64 = 0, length: Int32 = 0, amb: UInt8 = 4) {
        self.offset = offset
        self.length = length
        self.amb = amb
    }
}

/// Reference metadata: annotation + ambiguity information.
///
/// Corresponds to `bntseq_t` in bwa's `bntseq.h`.
public struct ReferenceMetadata: Sendable {
    /// Total length of the concatenated reference
    public var totalLength: Int64
    /// Number of reference sequences
    public var numSequences: Int32
    /// Seed for random (used for computing with Ns)
    public var seed: UInt32
    /// Per-sequence annotations
    public var annotations: [ReferenceAnnotation]
    /// Ambiguous regions
    public var ambiguities: [AmbiguityRegion]

    public init(totalLength: Int64 = 0, numSequences: Int32 = 0, seed: UInt32 = 11,
                annotations: [ReferenceAnnotation] = [], ambiguities: [AmbiguityRegion] = []) {
        self.totalLength = totalLength
        self.numSequences = numSequences
        self.seed = seed
        self.annotations = annotations
        self.ambiguities = ambiguities
    }

    /// Find the reference sequence ID for a given position in the concatenated reference.
    /// Uses binary search on annotation offsets.
    public func sequenceID(for position: Int64) -> Int32 {
        var lo = 0
        var hi = Int(numSequences) - 1
        while lo < hi {
            let mid = (lo + hi + 1) / 2
            if annotations[mid].offset <= position {
                lo = mid
            } else {
                hi = mid - 1
            }
        }
        return Int32(lo)
    }

    /// Convert a position in concatenated reference space to (rid, local_pos).
    public func decodePosition(_ pos: Int64) -> (rid: Int32, localPos: Int64) {
        let rid = sequenceID(for: pos)
        let localPos = pos - annotations[Int(rid)].offset
        return (rid, localPos)
    }
}
