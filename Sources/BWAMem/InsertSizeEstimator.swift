import BWACore

/// Running estimate of insert size distribution for paired-end reads.
/// Reimplements `mem_pestat()` from bwa-mem2's bwamem_pair.cpp.
public struct InsertSizeEstimator: Sendable {
    public private(set) var mean: Double = 0
    public private(set) var stddev: Double = 0
    public private(set) var low: Double = 0
    public private(set) var high: Double = 0
    public private(set) var count: Int = 0

    private var sum: Double = 0
    private var sumSq: Double = 0

    public init() {}

    /// Add an observed insert size.
    public mutating func add(_ insertSize: Int64) {
        let x = Double(abs(insertSize))
        count += 1
        sum += x
        sumSq += x * x

        mean = sum / Double(count)
        if count > 1 {
            let variance = (sumSq - sum * sum / Double(count)) / Double(count - 1)
            stddev = variance > 0 ? variance.squareRoot() : 0
        }

        // Set bounds at +/- 4 stddev
        low = mean - 4 * stddev
        high = mean + 4 * stddev
    }

    /// Check if an insert size is within expected range.
    public func isProperPair(_ insertSize: Int64) -> Bool {
        guard count >= 25 else { return true }  // Not enough data
        let x = Double(abs(insertSize))
        return x >= low && x <= high
    }
}
