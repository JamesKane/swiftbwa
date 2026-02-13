import Foundation
import BWACore
import CHelpers

/// Reads FASTA (plain or gzipped) and produces .pac, .ann, .amb files
/// matching bwa-mem2's format exactly.
struct FASTAEncoder: Sendable {

    struct EncodeResult {
        let pacLen: Int64
        let metadata: ReferenceMetadata
    }

    /// Encode a FASTA file into .pac, .ann, .amb files.
    /// Returns the forward-strand pac length and metadata.
    static func encode(fasta: String, prefix: String, verbosity: Int) throws -> EncodeResult {
        // Open FASTA (supports gzip via zlib)
        let gz = gzopen(fasta, "r")
        guard gz != nil else {
            throw BWAError.fileIOError("Cannot open FASTA: \(fasta)")
        }
        defer { gzclose(gz) }

        // State for PAC encoding
        var pacBuf = [UInt8]()  // pac bytes being built
        var pacLen: Int64 = 0   // total bases encoded (forward strand)

        // Annotation/ambiguity tracking
        var annotations = [ReferenceAnnotation]()
        var ambiguities = [AmbiguityRegion]()

        // Current sequence state
        var seqName = ""
        var seqComment = ""
        var seqOffset: Int64 = 0
        var seqLen: Int64 = 0
        var seqNAmb: Int32 = 0

        // Ambiguity run tracking
        var ambStart: Int64 = -1
        var ambLen: Int64 = 0
        var ambChar: UInt8 = 0  // ASCII of the ambiguity character

        // N→random base: replicate bwa-mem2's srand48(11) + lrand48()&3
        c_srand48(11)

        // Finalize current ambiguity run
        func flushAmbRun() {
            if ambStart >= 0 && ambLen > 0 {
                ambiguities.append(AmbiguityRegion(
                    offset: ambStart, length: Int32(ambLen), amb: ambChar
                ))
                seqNAmb += 1
                ambStart = -1
                ambLen = 0
            }
        }

        // Finalize current sequence
        func flushSequence() {
            flushAmbRun()
            if seqLen > 0 {
                annotations.append(ReferenceAnnotation(
                    offset: seqOffset, length: Int32(seqLen), nAmb: seqNAmb,
                    gi: 0, name: seqName, anno: seqComment, isAlt: false
                ))
            }
        }

        // Pack one base into pacBuf
        func packBase(_ base: UInt8) {
            let idx = Int(pacLen >> 2)
            let shift = ((~Int(pacLen)) & 3) << 1
            if idx >= pacBuf.count {
                pacBuf.append(0)
            }
            pacBuf[idx] |= base << shift
            pacLen += 1
            seqLen += 1
        }

        // Read line by line
        let bufSize = 65536
        let lineBuf = UnsafeMutablePointer<CChar>.allocate(capacity: bufSize)
        defer { lineBuf.deallocate() }

        while gzgets(gz, lineBuf, Int32(bufSize)) != nil {
            // Trim trailing newline/whitespace
            var len = strlen(lineBuf)
            while len > 0 && (lineBuf[len - 1] == 0x0A || lineBuf[len - 1] == 0x0D || lineBuf[len - 1] == 0x20) {
                len -= 1
            }
            guard len > 0 else { continue }

            if lineBuf[0] == 0x3E {  // '>'
                // Header line
                flushSequence()

                // Parse name and comment from header
                let header: String
                if len > 1 {
                    header = String(
                        data: Data(bytes: lineBuf + 1, count: len - 1),
                        encoding: .ascii
                    ) ?? ""
                } else {
                    header = ""
                }

                let trimmed = header.trimmingCharacters(in: .whitespaces)
                if let spaceIdx = trimmed.firstIndex(where: { $0 == " " || $0 == "\t" }) {
                    seqName = String(trimmed[..<spaceIdx])
                    seqComment = String(trimmed[trimmed.index(after: spaceIdx)...])
                } else {
                    seqName = trimmed
                    seqComment = ""
                }

                seqOffset = pacLen
                seqLen = 0
                seqNAmb = 0
                ambStart = -1
                ambLen = 0

                if verbosity >= 3 {
                    fputs("[M::index] Processing sequence: \(seqName)\n", stderr)
                }
            } else {
                // Sequence line
                for i in 0..<len {
                    let ch = lineBuf[i]
                    let base: UInt8
                    let isAmb: Bool

                    switch ch {
                    case 0x41, 0x61: base = 0; isAmb = false  // A/a
                    case 0x43, 0x63: base = 1; isAmb = false  // C/c
                    case 0x47, 0x67: base = 2; isAmb = false  // G/g
                    case 0x54, 0x74: base = 3; isAmb = false  // T/t
                    default:
                        // Ambiguous base → random via lrand48
                        base = UInt8(c_lrand48() & 3)
                        isAmb = true
                    }

                    if isAmb {
                        let upperCh = UInt8(ch >= 0x61 && ch <= 0x7A ? ch - 32 : ch)
                        if ambStart >= 0 && upperCh == ambChar {
                            ambLen += 1
                        } else {
                            flushAmbRun()
                            ambStart = pacLen
                            ambLen = 1
                            ambChar = upperCh
                        }
                    } else {
                        flushAmbRun()
                    }

                    packBase(base)
                }
            }
        }

        // Finalize last sequence
        flushSequence()

        // Write .pac file
        // Format: pac_bytes + zero_byte_if_aligned + remainder_byte
        let pacPath = prefix + ".pac"
        guard let pacFile = fopen(pacPath, "wb") else {
            throw BWAError.fileIOError("Cannot create: \(pacPath)")
        }
        defer { fclose(pacFile) }

        if !pacBuf.isEmpty {
            _ = pacBuf.withUnsafeBufferPointer { buf in
                fwrite(buf.baseAddress!, 1, buf.count, pacFile)
            }
        }
        // If pacLen%4 == 0, write an extra zero byte before the remainder
        let remainder = UInt8(pacLen % 4)
        if remainder == 0 {
            var zero: UInt8 = 0
            fwrite(&zero, 1, 1, pacFile)
        }
        var rem = remainder
        fwrite(&rem, 1, 1, pacFile)

        // Write .ann file
        let annPath = prefix + ".ann"
        guard let annFile = fopen(annPath, "w") else {
            throw BWAError.fileIOError("Cannot create: \(annPath)")
        }
        defer { fclose(annFile) }

        fputs("\(pacLen) \(annotations.count) 11\n", annFile)
        for ann in annotations {
            if ann.anno.isEmpty {
                fputs("\(ann.gi) \(ann.name)\n", annFile)
            } else {
                fputs("\(ann.gi) \(ann.name) \(ann.anno)\n", annFile)
            }
            fputs("\(ann.offset) \(ann.length) \(ann.nAmb)\n", annFile)
        }

        // Write .amb file
        let ambPath = prefix + ".amb"
        guard let ambFile = fopen(ambPath, "w") else {
            throw BWAError.fileIOError("Cannot create: \(ambPath)")
        }
        defer { fclose(ambFile) }

        fputs("\(pacLen) \(annotations.count) \(ambiguities.count)\n", ambFile)
        for amb in ambiguities {
            fputs("\(amb.offset) \(amb.length) \(UnicodeScalar(amb.amb))\n", ambFile)
        }

        let metadata = ReferenceMetadata(
            totalLength: pacLen,
            numSequences: Int32(annotations.count),
            seed: 11,
            annotations: annotations,
            ambiguities: ambiguities
        )

        return EncodeResult(pacLen: pacLen, metadata: metadata)
    }
}
