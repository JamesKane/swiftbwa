import BWACore
import Htslib

/// Converts alignment regions to BAM records using swift-htslib.
public struct SAMOutputBuilder: Sendable {

    /// Build a SAM header from reference metadata.
    public static func buildHeader(metadata: ReferenceMetadata) throws -> SAMHeader {
        let header = try SAMHeader()

        // @HD line
        try header.addLine(type: "HD", keyValues: [("VN", "1.6"), ("SO", "unsorted")])

        // @SQ lines for each reference sequence
        for ann in metadata.annotations {
            try header.addLine(type: "SQ", keyValues: [
                ("SN", ann.name),
                ("LN", String(ann.length))
            ])
        }

        // @PG line
        try header.addLine(type: "PG", keyValues: [
            ("ID", "swiftbwa"),
            ("PN", "swiftbwa"),
            ("VN", "0.1.0")
        ])

        return header
    }

    /// Convert an alignment region to a BAM record.
    ///
    /// The returned `BAMRecord` is move-only (`~Copyable`). Callers must consume it
    /// (e.g., by writing it to a file) rather than copying.
    public static func buildRecord(
        read: ReadSequence,
        region: MemAlnReg,
        allRegions: [MemAlnReg],
        metadata: ReferenceMetadata,
        scoring: ScoringParameters,
        cigar: [UInt32],
        nm: Int32 = 0,
        isPrimary: Bool,
        adjustedPos: Int64? = nil,
        mateRecord: (tid: Int32, pos: Int64, isReverse: Bool)? = nil
    ) throws -> BAMRecord {
        let mapq = MappingQuality.compute(
            region: region,
            allRegions: allRegions,
            scoring: scoring,
            readLength: Int32(read.length)
        )

        // Build flags
        var flag: AlignmentFlag = []
        if mateRecord != nil {
            flag.insert(.paired)
        }
        if region.rb < metadata.totalLength && region.re > 0 {
            // Check if reverse strand: position >= genome length means reverse
            let genomeLen = metadata.totalLength
            if region.rb >= genomeLen {
                flag.insert(.reverse)
            }
        }
        if !isPrimary {
            flag.insert(.secondary)
        }
        if let mate = mateRecord {
            if mate.isReverse {
                flag.insert(.mateReverse)
            }
        }

        // Compute position in reference coordinates
        // Use adjustedPos if provided (accounts for leading deletion squeeze)
        let refPos = adjustedPos ?? region.rb
        let (rid, localPos) = metadata.decodePosition(refPos)

        // Convert sequence to ASCII string
        let seqStr = String(read.bases.map { b -> Character in
            switch b {
            case 0: return "A"
            case 1: return "C"
            case 2: return "G"
            case 3: return "T"
            default: return "N"
            }
        })

        // Convert qualities to ASCII string (Phred+33)
        let qualStr = String(read.qualities.map { Character(UnicodeScalar($0 + 33)) })

        let mtid: Int32 = mateRecord?.tid ?? -1
        let mpos: Int64 = mateRecord?.pos ?? -1
        let isize: Int64 = 0

        var record = try BAMRecord()
        try record.set(
            qname: read.name,
            flag: flag.rawValue,
            tid: rid,
            pos: localPos,
            mapq: mapq,
            cigar: cigar,
            mtid: mtid,
            mpos: mpos,
            isize: isize,
            seq: seqStr,
            qual: qualStr
        )

        // Add auxiliary tags
        let aux = record.mutableAuxiliaryData
        try aux.updateInt(tag: "AS", value: Int64(region.score))
        try aux.updateInt(tag: "XS", value: Int64(region.sub))
        try aux.updateInt(tag: "NM", value: Int64(nm))

        return record
    }

    /// Build an unmapped BAM record.
    ///
    /// The returned `BAMRecord` is move-only (`~Copyable`). Callers must consume it
    /// (e.g., by writing it to a file) rather than copying.
    public static func buildUnmappedRecord(read: ReadSequence) throws -> BAMRecord {
        let seqStr = String(read.bases.map { b -> Character in
            switch b {
            case 0: return "A"
            case 1: return "C"
            case 2: return "G"
            case 3: return "T"
            default: return "N"
            }
        })
        let qualStr = String(read.qualities.map { Character(UnicodeScalar($0 + 33)) })

        var record = try BAMRecord()
        try record.set(
            qname: read.name,
            flag: AlignmentFlag.unmapped.rawValue,
            tid: -1,
            pos: -1,
            mapq: 0,
            cigar: [],
            mtid: -1,
            mpos: -1,
            isize: 0,
            seq: seqStr,
            qual: qualStr
        )

        return record
    }
}
