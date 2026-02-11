import BWACore

/// User-facing configuration for BWA-MEM alignment.
public struct BWAMemOptions: Sendable {
    public var scoring: ScoringParameters
    public var outputMode: OutputMode
    public var isPairedEnd: Bool

    public enum OutputMode: Sendable {
        case sam
        case bam
    }

    public init() {
        self.scoring = ScoringParameters()
        self.outputMode = .sam
        self.isPairedEnd = false
    }
}
