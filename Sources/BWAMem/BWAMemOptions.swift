import BWACore

/// Manual insert size override values (-I flag).
public struct InsertSizeOverride: Sendable {
    public var mean: Double
    public var stddev: Double
    public var max: Double
    public var min: Double

    public init(mean: Double, stddev: Double, max: Double? = nil, min: Double? = nil) {
        self.mean = mean
        self.stddev = stddev
        self.max = max ?? (mean + 3.5 * stddev)
        self.min = min ?? Swift.max(mean - 7.0 * stddev, 1.0)
    }
}

/// User-facing configuration for BWA-MEM alignment.
public struct BWAMemOptions: Sendable {
    public var scoring: ScoringParameters
    public var isPairedEnd: Bool
    public var readGroupLine: String?
    /// Extra header lines to insert (-H)
    public var headerLines: String?
    /// Append FASTQ comment to SAM output (-C)
    public var appendComment: Bool = false
    /// Skip ALT contig loading (-j)
    public var ignoreAlt: Bool = false
    /// Output reference header in XR tag (-V)
    public var outputRefHeader: Bool = false
    /// Manual insert size override (-I)
    public var manualInsertSize: InsertSizeOverride? = nil
    /// Verbosity level (-v): 1=error, 2=warning, 3=info (default), 4+=debug
    public var verbosity: Int = 3
    /// Use Metal GPU acceleration for Smith-Waterman kernels
    public var useGPU: Bool = false
    /// Emit Z-prefix training tags for MAPQ model training
    public var emitTrainingTags: Bool = false

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
        self.isPairedEnd = false
        self.readGroupLine = nil
    }
}
