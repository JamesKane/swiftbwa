import Foundation
import BWACore

public struct FMIndexLoader: Sendable {

    /// Load a complete FM-Index from the given prefix path.
    /// Expects files: {prefix}.bwt.2bit.64, {prefix}.ann, {prefix}.amb, {prefix}.pac
    public static func load(from prefix: String, skipAlt: Bool = false) throws -> FMIndex {
        let (bwt, bwtMappedFile, refSeqLen) = try loadBWT(from: prefix)
        let sa = try loadSA(from: prefix, referenceSeqLen: refSeqLen, mappedFile: bwtMappedFile)
        let pac = try loadPac(from: prefix)
        var metadata = try loadMetadata(from: prefix)
        if !skipAlt {
            loadAlt(from: prefix, metadata: &metadata)
        }

        return FMIndex(bwt: bwt, suffixArray: sa, packedRef: pac, metadata: metadata)
    }

    // MARK: - BWT Loading (.bwt.2bit.64)

    /// Load BWT from .bwt.2bit.64 file via mmap (matches FMI_search.cpp:384-494)
    /// Returns the BWT, the shared MappedFile (reused by SA), and referenceSeqLen.
    static func loadBWT(from prefix: String) throws -> (bwt: BWT, mappedFile: MappedFile, referenceSeqLen: Int64) {
        let path = prefix + ".bwt.2bit.64"
        let mf = try MappedFile(path: path)
        let base = mf.pointer

        // offset 0: referenceSeqLen (8 bytes)
        let referenceSeqLen = base.load(as: Int64.self)
        guard referenceSeqLen > 0 else {
            throw BWAError.indexCorrupted("Invalid reference_seq_len: \(referenceSeqLen)")
        }

        // offset 8: count[5] (40 bytes) — add 1 to each (matching bwa-mem2 load_index lines 432-436)
        let countPtr = (base + 8).assumingMemoryBound(to: Int64.self)
        let count = (countPtr[0] + 1, countPtr[1] + 1, countPtr[2] + 1, countPtr[3] + 1, countPtr[4] + 1)

        // offset 48: checkpoint OCC array
        let cpOccCount = Int(referenceSeqLen >> CP_SHIFT) + 1
        let cpStart = base + 48
        let cpBuffer = cpStart.assumingMemoryBound(to: CheckpointOCC.self)
        let checkpoints = UnsafeBufferPointer(start: cpBuffer, count: cpOccCount)

        // sentinel_index: last 8 bytes of file (may not be 8-byte aligned)
        let sentinelIndex = (base + mf.size - 8).loadUnaligned(as: Int64.self)

        let bwt = BWT(
            checkpoints: checkpoints,
            count: count,
            length: referenceSeqLen,
            sentinelIndex: sentinelIndex,
            mappedFile: mf
        )
        return (bwt, mf, referenceSeqLen)
    }

    // MARK: - SA Loading (.bwt.2bit.64)

    /// Load suffix array from the same .bwt.2bit.64 mmap'd file.
    /// Reuses the MappedFile from loadBWT — both BWT and SA hold strong references.
    static func loadSA(from prefix: String, referenceSeqLen: Int64, mappedFile mf: MappedFile) throws -> SuffixArray {
        let base = mf.pointer

        let cpOccCount = Int(referenceSeqLen >> CP_SHIFT) + 1
        let cpOccByteSize = cpOccCount * MemoryLayout<CheckpointOCC>.size
        let saOffset = 48 + cpOccByteSize

        let saCompX = 3
        let saCount = Int(referenceSeqLen >> saCompX) + 1

        let msStart = (base + saOffset).assumingMemoryBound(to: Int8.self)
        // ls_words may not be 4-byte aligned (follows variable-length ms_bytes)
        let lsRaw = base + saOffset + saCount

        return SuffixArray(
            msBytes: UnsafeBufferPointer(start: msStart, count: saCount),
            lsRawBase: lsRaw,
            count: saCount,
            compressionShift: saCompX,
            mappedFile: mf
        )
    }

    // MARK: - .pac Loading

    /// Load packed reference from .pac file via mmap.
    public static func loadPac(from prefix: String) throws -> PackedReference {
        let path = prefix + ".pac"
        let mf = try MappedFile(path: path)
        let base = mf.pointer

        guard mf.size >= 1 else {
            throw BWAError.indexCorrupted("Empty .pac file")
        }

        // Last byte stores remainder count
        let lastByte = (base + mf.size - 1).load(as: UInt8.self)
        let byteCount = mf.size - 1
        let pacLen: Int64
        if lastByte == 0 {
            pacLen = Int64(byteCount) * 4
        } else {
            pacLen = (Int64(byteCount) - 1) * 4 + Int64(lastByte)
        }

        let dataPtr = base.assumingMemoryBound(to: UInt8.self)
        return PackedReference(
            data: UnsafeBufferPointer(start: dataPtr, count: byteCount),
            length: pacLen,
            mappedFile: mf
        )
    }

    // MARK: - .alt Loading

    /// Load ALT contig annotations from .alt file.
    /// Lines starting with `@` are headers and skipped. Each other line's first
    /// field (before tab) is matched against annotation names to set `isAlt`.
    public static func loadAlt(from prefix: String, metadata: inout ReferenceMetadata) {
        let altPath = prefix + ".alt"
        guard let content = try? String(contentsOfFile: altPath, encoding: .utf8) else { return }

        // Build name -> index lookup
        var nameToIdx: [String: Int] = [:]
        for (i, ann) in metadata.annotations.enumerated() {
            nameToIdx[ann.name] = i
        }

        // Parse .alt file: each line's first field (before tab) is a sequence name
        for line in content.split(separator: "\n") {
            let name = String(line.prefix(while: { $0 != "\t" && $0 != "\r" }))
            guard !name.isEmpty && !name.hasPrefix("@") else { continue }
            if let idx = nameToIdx[name] {
                metadata.annotations[idx].isAlt = true
            }
        }
    }

    // MARK: - .ann/.amb Loading

    public static func loadMetadata(from prefix: String) throws -> ReferenceMetadata {
        var metadata = ReferenceMetadata()

        // Load .ann file
        let annPath = prefix + ".ann"
        let annContent: String
        do {
            annContent = try String(contentsOfFile: annPath, encoding: .utf8)
        } catch {
            throw BWAError.indexNotFound(annPath)
        }

        let annLines = annContent.split(separator: "\n", omittingEmptySubsequences: false)
        guard annLines.count >= 1 else {
            throw BWAError.indexCorrupted("Empty .ann file")
        }

        // First line: l_pac n_seqs seed
        let headerParts = annLines[0].split(separator: " ")
        guard headerParts.count >= 3 else {
            throw BWAError.indexCorrupted("Invalid .ann header")
        }
        metadata.totalLength = Int64(headerParts[0]) ?? 0
        metadata.numSequences = Int32(headerParts[1]) ?? 0
        metadata.seed = UInt32(headerParts[2]) ?? 11

        // Each sequence has 2 lines: (gi offset name anno) and (length nAmb)
        var lineIdx = 1
        for _ in 0..<metadata.numSequences {
            guard lineIdx + 1 < annLines.count else { break }

            var ann = ReferenceAnnotation()
            let parts1 = annLines[lineIdx].split(separator: " ", maxSplits: 2)
            if parts1.count >= 2 {
                ann.gi = UInt32(parts1[0]) ?? 0
                ann.name = String(parts1[1])
                ann.anno = parts1.count > 2 ? String(parts1[2]) : ""
            }
            lineIdx += 1

            let parts2 = annLines[lineIdx].split(separator: " ")
            if parts2.count >= 3 {
                ann.offset = Int64(parts2[0]) ?? 0
                ann.length = Int32(parts2[1]) ?? 0
                ann.nAmb = Int32(parts2[2]) ?? 0
            }
            lineIdx += 1

            metadata.annotations.append(ann)
        }

        // Load .amb file
        let ambPath = prefix + ".amb"
        if let ambContent = try? String(contentsOfFile: ambPath, encoding: .utf8) {
            let ambLines = ambContent.split(separator: "\n", omittingEmptySubsequences: false)
            if ambLines.count >= 1 {
                // First line: l_pac n_seqs n_holes
                let ambHeader = ambLines[0].split(separator: " ")
                let nHoles = ambHeader.count >= 3 ? (Int(ambHeader[2]) ?? 0) : 0

                var ambLineIdx = 1
                for _ in 0..<nHoles {
                    guard ambLineIdx < ambLines.count else { break }
                    let parts = ambLines[ambLineIdx].split(separator: " ")
                    if parts.count >= 3 {
                        let region = AmbiguityRegion(
                            offset: Int64(parts[0]) ?? 0,
                            length: Int32(parts[1]) ?? 0,
                            amb: UInt8(parts[2].first?.asciiValue ?? 78)  // 'N'
                        )
                        metadata.ambiguities.append(region)
                    }
                    ambLineIdx += 1
                }
            }
        }

        return metadata
    }
}
