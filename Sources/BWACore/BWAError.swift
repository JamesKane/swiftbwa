import Foundation

public enum BWAError: Error, Sendable {
    case indexNotFound(String)
    case indexCorrupted(String)
    case invalidFileFormat(String)
    case fileIOError(String)
    case invalidSequence(String)
    case alignmentFailed(String)
    case invalidParameter(String)
    case mmapFailed(String)
}

extension BWAError: LocalizedError {
    public var errorDescription: String? {
        switch self {
        case .indexNotFound(let msg): return "Index not found: \(msg)"
        case .indexCorrupted(let msg): return "Index corrupted: \(msg)"
        case .invalidFileFormat(let msg): return "Invalid file format: \(msg)"
        case .fileIOError(let msg): return "File I/O error: \(msg)"
        case .invalidSequence(let msg): return "Invalid sequence: \(msg)"
        case .alignmentFailed(let msg): return "Alignment failed: \(msg)"
        case .invalidParameter(let msg): return "Invalid parameter: \(msg)"
        case .mmapFailed(let msg): return "Memory mapping failed: \(msg)"
        }
    }
}
