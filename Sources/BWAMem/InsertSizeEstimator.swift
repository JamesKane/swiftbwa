import BWACore
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif canImport(Musl)
import Musl
#endif

/// Orientation of a read pair.
public enum PairOrientation: Int, Sendable, CaseIterable {
    case fr = 0  // Forward-Reverse (standard Illumina)
    case rf = 1  // Reverse-Forward
    case ff = 2  // Forward-Forward
    case rr = 3  // Reverse-Reverse
}

/// Per-orientation insert size statistics.
public struct OrientationStats: Sendable {
    public var low: Int64 = 0
    public var high: Int64 = 0
    public var mean: Double = 0
    public var stddev: Double = 0
    public var count: Int = 0
    public var failed: Bool = true

    public init() {}
}

/// Insert size distribution estimated from high-quality concordant pairs.
public struct InsertSizeDistribution: Sendable {
    public var stats: [OrientationStats]
    public var primaryOrientation: PairOrientation

    public init() {
        stats = Array(repeating: OrientationStats(), count: 4)
        primaryOrientation = .fr
    }
}

/// Quartile-based insert size estimation following bwa-mem2's `mem_pestat()`.
public struct InsertSizeEstimator: Sendable {

    /// Minimum number of samples per orientation to compute stats.
    public static let minSamples = 25

    /// Infer pair orientation and unsigned insert size from two alignment regions.
    ///
    /// Uses the BWT coordinate convention: positions >= genomeLength are on the
    /// reverse-complement strand.
    ///
    /// - Parameters:
    ///   - r1: Primary alignment region for read 1
    ///   - r2: Primary alignment region for read 2
    ///   - genomeLength: Half of the total BWT length (forward genome length)
    /// - Returns: Tuple of (orientation, unsigned insert size), or nil if on different chromosomes
    public static func inferOrientation(
        r1: MemAlnReg, r2: MemAlnReg, genomeLength: Int64
    ) -> (orientation: PairOrientation, insertSize: Int64)? {
        guard r1.rid == r2.rid else { return nil }

        let r1Reverse = r1.rb >= genomeLength
        let r2Reverse = r2.rb >= genomeLength

        // Convert to forward-strand positions
        let r1Start: Int64
        let r1End: Int64
        if r1Reverse {
            r1Start = 2 * genomeLength - r1.re
            r1End = 2 * genomeLength - r1.rb
        } else {
            r1Start = r1.rb
            r1End = r1.re
        }

        let r2Start: Int64
        let r2End: Int64
        if r2Reverse {
            r2Start = 2 * genomeLength - r2.re
            r2End = 2 * genomeLength - r2.rb
        } else {
            r2Start = r2.rb
            r2End = r2.re
        }

        // Insert size = distance between outer coordinates
        let outerLeft = min(r1Start, r2Start)
        let outerRight = max(r1End, r2End)
        let isize = outerRight - outerLeft

        // Determine orientation based on strand and relative position
        let orientation: PairOrientation
        if !r1Reverse && r2Reverse {
            // Read1 forward, read2 reverse
            orientation = r1Start <= r2Start ? .fr : .rf
        } else if r1Reverse && !r2Reverse {
            // Read1 reverse, read2 forward
            orientation = r2Start <= r1Start ? .fr : .rf
        } else if !r1Reverse && !r2Reverse {
            orientation = .ff
        } else {
            orientation = .rr
        }

        return (orientation, isize)
    }

    /// Estimate insert size distribution from a set of high-quality pairs.
    ///
    /// Implements bwa-mem2's `mem_pestat()` quartile-based approach:
    /// 1. Collect insert sizes per orientation from pairs with MAPQ >= 20
    /// 2. For each orientation with >= 25 samples:
    ///    - Compute Q25, Q75 and IQR
    ///    - Filter to [Q25 - 2*IQR, Q75 + 2*IQR]
    ///    - Compute mean and stddev on filtered set
    /// 3. Pick primary orientation as the one with most samples
    ///
    /// - Parameters:
    ///   - regions1: Alignment regions for each read1 (indexed by read pair index)
    ///   - regions2: Alignment regions for each read2 (indexed by read pair index)
    ///   - genomeLength: Forward genome length
    /// - Returns: Insert size distribution
    public static func estimate(
        regions1: [[MemAlnReg]],
        regions2: [[MemAlnReg]],
        genomeLength: Int64
    ) -> InsertSizeDistribution {
        var dist = InsertSizeDistribution()

        // Collect insert sizes per orientation
        var samples: [[Int64]] = [[], [], [], []]  // indexed by PairOrientation.rawValue

        let pairCount = min(regions1.count, regions2.count)
        for i in 0..<pairCount {
            // Find primary region for each end (secondary < 0)
            guard let r1 = regions1[i].first(where: { $0.secondary < 0 }),
                  let r2 = regions2[i].first(where: { $0.secondary < 0 }) else {
                continue
            }

            // Require reasonable mapping quality: score should be significantly
            // better than sub-optimal
            let scoring = ScoringParameters()
            let mapq1 = MappingQuality.compute(
                region: r1, allRegions: regions1[i],
                scoring: scoring, readLength: r1.qe
            )
            let mapq2 = MappingQuality.compute(
                region: r2, allRegions: regions2[i],
                scoring: scoring, readLength: r2.qe
            )
            guard mapq1 >= 20 && mapq2 >= 20 else { continue }

            guard let result = inferOrientation(r1: r1, r2: r2, genomeLength: genomeLength) else {
                continue
            }

            samples[result.orientation.rawValue].append(result.insertSize)
        }

        // Compute stats per orientation
        var maxCount = 0
        var primaryIdx = 0

        for ori in PairOrientation.allCases {
            let idx = ori.rawValue
            var sizes = samples[idx]
            guard sizes.count >= Self.minSamples else { continue }

            sizes.sort()

            // Compute Q25 and Q75
            let q25 = sizes[sizes.count / 4]
            let q75 = sizes[(sizes.count * 3) / 4]
            let iqr = q75 - q25

            // Set bounds
            let low = q25 - 2 * iqr
            let high = q75 + 2 * iqr

            // Filter to within bounds and recompute
            let filtered = sizes.filter { $0 >= low && $0 <= high }
            guard !filtered.isEmpty else { continue }

            let sum = filtered.reduce(0.0) { $0 + Double($1) }
            let mean = sum / Double(filtered.count)
            let variance = filtered.reduce(0.0) { $0 + (Double($1) - mean) * (Double($1) - mean) }
                / Double(max(1, filtered.count - 1))
            let stddev = sqrt(variance)

            dist.stats[idx].low = low
            dist.stats[idx].high = high
            dist.stats[idx].mean = mean
            dist.stats[idx].stddev = stddev
            dist.stats[idx].count = filtered.count
            dist.stats[idx].failed = false

            if filtered.count > maxCount {
                maxCount = filtered.count
                primaryIdx = idx
            }
        }

        dist.primaryOrientation = PairOrientation(rawValue: primaryIdx) ?? .fr
        return dist
    }
}
