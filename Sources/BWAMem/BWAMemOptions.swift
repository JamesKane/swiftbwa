import BWACore

/// User-facing configuration for BWA-MEM alignment.
public struct BWAMemOptions: Sendable {
    public var scoring: ScoringParameters
    public var outputMode: OutputMode
    public var isPairedEnd: Bool
    public var readGroupLine: String?
    /// Extra header lines to insert (-H)
    public var headerLines: String?
    /// Append FASTQ comment to SAM output (-C)
    public var appendComment: Bool = false
    /// Skip ALT contig loading (-j)
    public var ignoreAlt: Bool = false

    public enum OutputMode: Sendable {
        case sam
        case bam
    }

    /// Extract the ID field from the read group line.
    public var readGroupID: String? {
        guard let line = readGroupLine else { return nil }
        let fields = line.replacingOccurrences(of: "\\t", with: "\t")
                         .split(separator: "\t")
        for field in fields {
            if field.hasPrefix("ID:") {
                return String(field.dropFirst(3))
            }
        }
        return nil
    }

    public init() {
        self.scoring = ScoringParameters()
        self.outputMode = .sam
        self.isPairedEnd = false
        self.readGroupLine = nil
    }
}
