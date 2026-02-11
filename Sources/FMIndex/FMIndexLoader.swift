import Foundation
import BWACore

public struct FMIndexLoader: Sendable {

    /// Load a complete FM-Index from the given prefix path.
    /// Expects files: {prefix}.bwt.2bit.64, {prefix}.ann, {prefix}.amb, {prefix}.pac
    public static func load(from prefix: String) throws -> FMIndex {
        let bwtResult = try loadBWT(from: prefix)
        let sa = try loadSA(from: prefix, referenceSeqLen: bwtResult.length)
        let pac = try loadPac(from: prefix)
        let metadata = try loadMetadata(from: prefix)

        return FMIndex(bwt: bwtResult, suffixArray: sa, packedRef: pac, metadata: metadata)
    }

    // MARK: - BWT Loading (.bwt.2bit.64)

    /// Load BWT from .bwt.2bit.64 file (matches FMI_search.cpp:384-494)
    public static func loadBWT(from prefix: String) throws -> BWT {
        let path = prefix + ".bwt.2bit.64"

        guard let fh = FileHandle(forReadingAtPath: path) else {
            throw BWAError.indexNotFound(path)
        }
        defer { fh.closeFile() }

        // Read reference_seq_len (int64)
        guard let lenData = try fh.read(upToCount: 8), lenData.count == 8 else {
            throw BWAError.indexCorrupted("Cannot read reference_seq_len")
        }
        let referenceSeqLen = lenData.withUnsafeBytes { $0.load(as: Int64.self) }
        guard referenceSeqLen > 0 else {
            throw BWAError.indexCorrupted("Invalid reference_seq_len: \(referenceSeqLen)")
        }

        // Read count[5] (5 x int64)
        guard let countData = try fh.read(upToCount: 40), countData.count == 40 else {
            throw BWAError.indexCorrupted("Cannot read count array")
        }
        var rawCount = (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0))
        countData.withUnsafeBytes { buf in
            let ptr = buf.bindMemory(to: Int64.self)
            rawCount = (ptr[0], ptr[1], ptr[2], ptr[3], ptr[4])
        }
        // Add 1 to each count (matching bwa-mem2 load_index lines 432-436)
        let count = (rawCount.0 + 1, rawCount.1 + 1, rawCount.2 + 1, rawCount.3 + 1, rawCount.4 + 1)

        // Read checkpoint OCC array
        let cpOccCount = Int(referenceSeqLen >> CP_SHIFT) + 1
        let cpOccByteSize = cpOccCount * MemoryLayout<CheckpointOCC>.size
        guard let cpData = try fh.read(upToCount: cpOccByteSize), cpData.count == cpOccByteSize else {
            throw BWAError.indexCorrupted("Cannot read checkpoint OCC array")
        }

        // Copy checkpoint data into allocated buffer
        let cpBuffer = UnsafeMutablePointer<CheckpointOCC>.allocate(capacity: cpOccCount)
        cpData.withUnsafeBytes { raw in
            let src = raw.bindMemory(to: CheckpointOCC.self)
            cpBuffer.initialize(from: src.baseAddress!, count: cpOccCount)
        }

        // Read sentinel index (at the end of BWT section, after SA data)
        // Actually: the .bwt.2bit.64 layout is:
        //   reference_seq_len (8 bytes)
        //   count[5] (40 bytes)
        //   cp_occ[] (cpOccByteSize bytes)
        //   sa_ms_byte[] (saCount bytes)
        //   sa_ls_word[] (saCount * 4 bytes)
        //   sentinel_index (8 bytes)
        // We skip SA data to read sentinel index
        let saCompX = 3
        let saCount = Int(referenceSeqLen >> saCompX) + 1
        let saSkipBytes = saCount + saCount * 4
        guard let _ = try fh.read(upToCount: saSkipBytes) else {
            throw BWAError.indexCorrupted("Cannot skip SA data to read sentinel_index")
        }

        guard let sentData = try fh.read(upToCount: 8), sentData.count == 8 else {
            throw BWAError.indexCorrupted("Cannot read sentinel_index")
        }
        let sentinelIndex = sentData.withUnsafeBytes { $0.load(as: Int64.self) }

        return BWT(
            checkpoints: UnsafeBufferPointer(start: cpBuffer, count: cpOccCount),
            count: count,
            length: referenceSeqLen,
            sentinelIndex: sentinelIndex,
            ownedBase: UnsafeMutableRawPointer(cpBuffer)
        )
    }

    // MARK: - SA Loading (.bwt.2bit.64)

    /// Load suffix array from .bwt.2bit.64 file
    public static func loadSA(from prefix: String, referenceSeqLen: Int64) throws -> SuffixArray {
        let path = prefix + ".bwt.2bit.64"
        guard let fh = FileHandle(forReadingAtPath: path) else {
            throw BWAError.indexNotFound(path)
        }
        defer { fh.closeFile() }

        // Skip: reference_seq_len (8) + count[5] (40) + checkpoints
        let cpOccCount = Int(referenceSeqLen >> CP_SHIFT) + 1
        let cpOccByteSize = cpOccCount * MemoryLayout<CheckpointOCC>.size
        let skipBytes = UInt64(8 + 40 + cpOccByteSize)
        fh.seek(toFileOffset: skipBytes)

        let saCompX = 3
        let saCount = Int(referenceSeqLen >> saCompX) + 1

        let msBytesSize = saCount
        guard let msData = try fh.read(upToCount: msBytesSize), msData.count == msBytesSize else {
            throw BWAError.indexCorrupted("Cannot read SA ms_byte")
        }

        let lsWordsSize = saCount * 4
        guard let lsData = try fh.read(upToCount: lsWordsSize), lsData.count == lsWordsSize else {
            throw BWAError.indexCorrupted("Cannot read SA ls_word")
        }

        let msBuffer = UnsafeMutablePointer<Int8>.allocate(capacity: saCount)
        msData.withUnsafeBytes { raw in
            raw.baseAddress!.withMemoryRebound(to: Int8.self, capacity: saCount) { src in
                msBuffer.initialize(from: src, count: saCount)
            }
        }

        let lsBuffer = UnsafeMutablePointer<UInt32>.allocate(capacity: saCount)
        lsData.withUnsafeBytes { raw in
            raw.baseAddress!.withMemoryRebound(to: UInt32.self, capacity: saCount) { src in
                lsBuffer.initialize(from: src, count: saCount)
            }
        }

        return SuffixArray(
            msBytes: UnsafeBufferPointer(start: msBuffer, count: saCount),
            lsWords: UnsafeBufferPointer(start: lsBuffer, count: saCount),
            count: saCount,
            compressionShift: saCompX,
            ownedMSBase: UnsafeMutableRawPointer(msBuffer),
            ownedLSBase: UnsafeMutableRawPointer(lsBuffer)
        )
    }

    // MARK: - .pac Loading

    public static func loadPac(from prefix: String) throws -> PackedReference {
        let path = prefix + ".pac"
        let url = URL(fileURLWithPath: path)
        let data: Data
        do {
            data = try Data(contentsOf: url)
        } catch {
            throw BWAError.indexNotFound(path)
        }
        guard data.count >= 1 else {
            throw BWAError.indexCorrupted("Empty .pac file")
        }

        // Last byte stores remainder count
        let lastByte = data[data.count - 1]
        let pacLen: Int64
        if lastByte == 0 {
            pacLen = Int64(data.count - 1) * 4
        } else {
            pacLen = (Int64(data.count - 1) - 1) * 4 + Int64(lastByte)
        }

        let byteCount = data.count - 1
        let buffer = UnsafeMutablePointer<UInt8>.allocate(capacity: byteCount)
        data.withUnsafeBytes { raw in
            raw.baseAddress!.withMemoryRebound(to: UInt8.self, capacity: byteCount) { src in
                buffer.initialize(from: src, count: byteCount)
            }
        }

        return PackedReference(
            data: UnsafeBufferPointer(start: buffer, count: byteCount),
            length: pacLen,
            ownedBase: UnsafeMutableRawPointer(buffer)
        )
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
