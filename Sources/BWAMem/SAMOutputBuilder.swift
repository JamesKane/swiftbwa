import Alignment
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
    public static func buildHeader(
        metadata: ReferenceMetadata,
        readGroupLine: String? = nil,
        headerLines: String? = nil
    ) throws -> SAMHeader {
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

        // @RG line (verbatim from user, with literal \t expanded)
        if let rgLine = readGroupLine {
            let expanded = rgLine.replacingOccurrences(of: "\\t", with: "\t")
            try header.addLines(expanded)
        }

        // Custom header lines (-H)
        if let hLines = headerLines {
            let expanded = hLines.replacingOccurrences(of: "\\t", with: "\t")
            try header.addLines(expanded)
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
        isSupplementary: Bool = false,
        mapqOverride: UInt8? = nil,
        adjustedPos: Int64? = nil,
        pairedEnd: PairedEndInfo? = nil,
        saTag: String? = nil,
        xaTag: String? = nil,
        readGroupID: String? = nil,
        outputRefHeader: Bool = false,
        appendComment: Bool = false
    ) throws -> BAMRecord {
        let mapq = mapqOverride ?? MappingQuality.compute(
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

        if region.rb >= metadata.totalLength {
            flag.insert(.reverse)
        }
        if isSupplementary {
            flag.insert(.supplementary)
        } else if !isPrimary {
            flag.insert(.secondary)
        }

        // Compute position in reference coordinates
        let refPos = adjustedPos ?? region.rb
        let (rid, localPos) = metadata.decodePosition(refPos)

        // Determine CIGAR, SEQ, QUAL — apply hard-clip for supplementary if not soft-clip mode
        var finalCigar = cigar
        var seqStr: String
        var qualStr: String

        // SAM spec: when FLAG 0x10, SEQ = reverse complement, QUAL = reversed
        let isReverse = flag.contains(.reverse)
        let fullSeq: [Character]
        let fullQual: [Character]
        if isReverse {
            fullSeq = read.bases.reversed().map { b -> Character in
                switch b {
                case 0: return "T"  // complement of A
                case 1: return "G"  // complement of C
                case 2: return "C"  // complement of G
                case 3: return "A"  // complement of T
                default: return "N"
                }
            }
            fullQual = read.qualities.reversed().map { Character(UnicodeScalar($0 + 33)) }
        } else {
            fullSeq = read.bases.map { b -> Character in
                switch b {
                case 0: return "A"
                case 1: return "C"
                case 2: return "G"
                case 3: return "T"
                default: return "N"
                }
            }
            fullQual = read.qualities.map { Character(UnicodeScalar($0 + 33)) }
        }

        if isSupplementary && (scoring.flag & ScoringParameters.flagSoftClip) == 0 {
            let (hardCigar, trimLeft, trimRight) = convertToHardClip(cigar: cigar)
            finalCigar = hardCigar
            let seqEnd = fullSeq.count - trimRight
            let seqStart = trimLeft
            if seqStart < seqEnd {
                seqStr = String(fullSeq[seqStart..<seqEnd])
                qualStr = String(fullQual[seqStart..<seqEnd])
            } else {
                seqStr = String(fullSeq)
                qualStr = String(fullQual)
            }
        } else {
            seqStr = String(fullSeq)
            qualStr = String(fullQual)
        }

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
            cigar: finalCigar,
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
        if region.altSc > 0 {
            let paValue = Float(region.score) / Float(region.altSc)
            try aux.updateFloat(tag: "pa", value: paValue)
        }
        if let sa = saTag {
            try aux.updateString(tag: "SA", value: sa)
        }
        if let xa = xaTag {
            try aux.updateString(tag: "XA", value: xa)
        }
        if let rg = readGroupID {
            try aux.updateString(tag: "RG", value: rg)
        }
        if outputRefHeader && rid >= 0 && rid < metadata.annotations.count {
            let anno = metadata.annotations[Int(rid)].anno
            if !anno.isEmpty {
                try aux.updateString(tag: "XR", value: anno)
            }
        }
        if appendComment && !read.comment.isEmpty {
            try aux.updateString(tag: "CO", value: read.comment)
        }

        return record
    }

    // MARK: - Hard-clip conversion

    /// Convert leading/trailing soft-clips to hard-clips.
    /// Returns the new CIGAR and trim lengths so the caller can substring SEQ/QUAL.
    public static func convertToHardClip(cigar: [UInt32]) -> (cigar: [UInt32], trimLeft: Int, trimRight: Int) {
        guard !cigar.isEmpty else { return (cigar, 0, 0) }
        var result = cigar
        var trimLeft = 0
        var trimRight = 0

        // Check leading op
        let firstOp = result[0] & 0xF
        if firstOp == CIGAROp.softClip.rawValue {
            let len = result[0] >> 4
            trimLeft = Int(len)
            result[0] = len << 4 | CIGAROp.hardClip.rawValue
        }

        // Check trailing op
        if result.count > 1 {
            let lastIdx = result.count - 1
            let lastOp = result[lastIdx] & 0xF
            if lastOp == CIGAROp.softClip.rawValue {
                let len = result[lastIdx] >> 4
                trimRight = Int(len)
                result[lastIdx] = len << 4 | CIGAROp.hardClip.rawValue
            }
        }

        return (result, trimLeft, trimRight)
    }

    // MARK: - SA tag builder

    /// Build an SA auxiliary tag string from alignment segments.
    /// Format per segment: `rname,pos,strand,CIGAR,mapQ,NM;`
    /// The SA CIGAR always uses soft-clips (not hard-clips), per bwa-mem2 convention.
    public static func buildSATag(
        segments: [(rname: String, pos: Int64, isReverse: Bool,
                    cigarString: String, mapq: UInt8, nm: Int32)],
        excludeIndex: Int
    ) -> String {
        var parts: [String] = []
        for (i, seg) in segments.enumerated() {
            if i == excludeIndex { continue }
            let strand = seg.isReverse ? "-" : "+"
            // SAM positions are 1-based
            parts.append("\(seg.rname),\(seg.pos + 1),\(strand),\(seg.cigarString),\(seg.mapq),\(seg.nm);")
        }
        return parts.joined()
    }

    // MARK: - XA tag builder

    /// Build an XA auxiliary tag string from secondary alignment info.
    /// Format per hit: `chr,{+|-}pos,CIGAR,NM;`
    /// Returns nil if count exceeds `maxHits` (bwa-mem2 behavior: too many hits → no XA tag).
    public static func buildXATag(
        secondaries: [(rname: String, pos: Int64, isReverse: Bool,
                       cigarString: String, nm: Int32)],
        maxHits: Int = 5
    ) -> String? {
        guard secondaries.count <= maxHits else { return nil }
        var parts: [String] = []
        for sec in secondaries {
            let strand = sec.isReverse ? "-" : "+"
            // SAM positions are 1-based
            parts.append("\(sec.rname),\(strand)\(sec.pos + 1),\(sec.cigarString),\(sec.nm);")
        }
        return parts.joined()
    }

    /// Build an unmapped BAM record.
    ///
    /// The returned `BAMRecord` is move-only (`~Copyable`). Callers must consume it
    /// (e.g., by writing it to a file) rather than copying.
    public static func buildUnmappedRecord(
        read: ReadSequence,
        pairedEnd: PairedEndInfo? = nil,
        readGroupID: String? = nil,
        appendComment: Bool = false
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

        if let rg = readGroupID {
            let aux = record.mutableAuxiliaryData
            try aux.updateString(tag: "RG", value: rg)
        }
        if appendComment && !read.comment.isEmpty {
            let aux = record.mutableAuxiliaryData
            try aux.updateString(tag: "CO", value: read.comment)
        }

        return record
    }
}
