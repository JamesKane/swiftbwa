import BWACore
import Htslib

/// Paired-end information for building SAM records.
public struct PairedEndInfo: Sendable {
    public var isRead1: Bool
    public var isProperPair: Bool
    public var mateTid: Int32
    public var matePos: Int64
    public var mateIsReverse: Bool
    public var mateIsUnmapped: Bool
    public var tlen: Int64
    public var mateCigarString: String?

    public init(
        isRead1: Bool, isProperPair: Bool,
        mateTid: Int32, matePos: Int64,
        mateIsReverse: Bool, mateIsUnmapped: Bool,
        tlen: Int64, mateCigarString: String? = nil
    ) {
        self.isRead1 = isRead1
        self.isProperPair = isProperPair
        self.mateTid = mateTid
        self.matePos = matePos
        self.mateIsReverse = mateIsReverse
        self.mateIsUnmapped = mateIsUnmapped
        self.tlen = tlen
        self.mateCigarString = mateCigarString
    }
}

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
        md: String? = nil,
        isPrimary: Bool,
        adjustedPos: Int64? = nil,
        pairedEnd: PairedEndInfo? = nil
    ) throws -> BAMRecord {
        let mapq = MappingQuality.compute(
            region: region,
            allRegions: allRegions,
            scoring: scoring,
            readLength: Int32(read.length)
        )

        // Build flags
        var flag: AlignmentFlag = []

        // Paired-end flags
        if let pe = pairedEnd {
            flag.insert(.paired)
            if pe.isRead1 {
                flag.insert(.read1)
            } else {
                flag.insert(.read2)
            }
            if pe.isProperPair {
                flag.insert(.properPair)
            }
            if pe.mateIsReverse {
                flag.insert(.mateReverse)
            }
            if pe.mateIsUnmapped {
                flag.insert(.mateUnmapped)
            }
        }

        if region.rb < metadata.totalLength && region.re > 0 {
            let genomeLen = metadata.totalLength
            if region.rb >= genomeLen {
                flag.insert(.reverse)
            }
        }
        if !isPrimary {
            flag.insert(.secondary)
        }

        // Compute position in reference coordinates
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

        let mtid: Int32
        let mpos: Int64
        let isize: Int64

        if let pe = pairedEnd {
            mtid = pe.mateIsUnmapped ? -1 : pe.mateTid
            mpos = pe.mateIsUnmapped ? -1 : pe.matePos
            isize = pe.tlen
        } else {
            mtid = -1
            mpos = -1
            isize = 0
        }

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
        if let md = md {
            try aux.updateString(tag: "MD", value: md)
        }
        if let pe = pairedEnd, let mc = pe.mateCigarString {
            try aux.updateString(tag: "MC", value: mc)
        }

        return record
    }

    /// Build an unmapped BAM record.
    ///
    /// The returned `BAMRecord` is move-only (`~Copyable`). Callers must consume it
    /// (e.g., by writing it to a file) rather than copying.
    public static func buildUnmappedRecord(
        read: ReadSequence,
        pairedEnd: PairedEndInfo? = nil
    ) throws -> BAMRecord {
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

        var flag: AlignmentFlag = [.unmapped]
        let mtid: Int32
        let mpos: Int64
        let isize: Int64 = 0

        if let pe = pairedEnd {
            flag.insert(.paired)
            if pe.isRead1 {
                flag.insert(.read1)
            } else {
                flag.insert(.read2)
            }
            if pe.mateIsReverse {
                flag.insert(.mateReverse)
            }
            if pe.mateIsUnmapped {
                flag.insert(.mateUnmapped)
            }
            mtid = pe.mateIsUnmapped ? -1 : pe.mateTid
            mpos = pe.mateIsUnmapped ? -1 : pe.matePos
        } else {
            mtid = -1
            mpos = -1
        }

        var record = try BAMRecord()
        try record.set(
            qname: read.name,
            flag: flag.rawValue,
            tid: -1,
            pos: -1,
            mapq: 0,
            cigar: [],
            mtid: mtid,
            mpos: mpos,
            isize: isize,
            seq: seqStr,
            qual: qualStr
        )

        return record
    }
}
