import BWACore
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#elseif canImport(Musl)
import Musl
#endif

/// Orientation of a read pair.
/// Raw values match bwa-mem2's direction encoding used in mem_matesw/mem_infer_dir:
/// bit 0 = same/diff strand, XOR 3 if mate upstream → FF=0, FR=1, RF=2, RR=3.
public enum PairOrientation: Int, Sendable, CaseIterable {
    case ff = 0  // Forward-Forward
    case fr = 1  // Forward-Reverse (standard Illumina)
    case rf = 2  // Reverse-Forward
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
    /// Lower proper-pair bound for rescue: Q25 - 3*IQR, capped at mean - 4*stddev
    public var properLow: Int64 = 0
    /// Upper proper-pair bound for rescue: Q75 + 3*IQR, capped at mean + 4*stddev
    public var properHigh: Int64 = 0

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
    /// bwa-mem2: MIN_DIR_CNT = 10
    public static let minSamples = 10

    /// Minimum ratio of an orientation's count to the max orientation's count.
    /// Orientations below this ratio are marked as failed.
    /// bwa-mem2: MIN_DIR_RATIO = 0.05
    public static let minDirRatio = 0.05

    /// Build an insert size distribution from manual override values (-I flag).
    ///
    /// Constructs FR-orientation stats from user-provided values, matching bwa-mem2's
    /// behavior for `-I mean,stddev[,max[,min]]`.
    public static func buildManualDistribution(override ov: InsertSizeOverride) -> InsertSizeDistribution {
        var dist = InsertSizeDistribution()
        var stats = OrientationStats()
        stats.failed = false
        stats.mean = ov.mean
        stats.stddev = ov.stddev
        stats.low = Int64(ov.min)
        stats.high = Int64(ov.max)
        stats.count = 1000  // synthetic count so it doesn't look failed
        stats.properLow = Swift.max(1, Int64(ov.mean - 4.0 * ov.stddev + 0.5))
        stats.properHigh = Int64(ov.mean + 4.0 * ov.stddev + 0.5)
        dist.stats[PairOrientation.fr.rawValue] = stats
        dist.primaryOrientation = .fr
        return dist
    }

    /// Infer pair orientation and unsigned insert size from two alignment regions.
    ///
    /// Uses bwa-mem2's `mem_infer_dir` algorithm with start-to-start distance metric
    /// in BWT coordinate space. Positions >= genomeLength are on the reverse strand.
    ///
    /// The bwa direction encoding is:
    /// - bit 0: 0 = same strand, 1 = different strand
    /// - XOR with 3 if p2 <= b1 (i.e., mate is upstream)
    /// - Result: FF=0, FR=1, RF=2, RR=3
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

        // Mirror r2 onto r1's strand if they're on different strands
        let sameStrand = (r1.rb >= genomeLength) == (r2.rb >= genomeLength)
        let p2 = sameStrand ? r2.rb : (2 * genomeLength - 1 - r2.rb)

        // Start-to-start distance
        let dist = abs(p2 - r1.rb)

        // bwa direction encoding: same-strand→0, diff→1, XOR with 3 if p2<=b1
        let bwaDir = (sameStrand ? 0 : 1) ^ (p2 > r1.rb ? 0 : 3)

        // bwa dir encoding matches PairOrientation rawValues: FF=0, FR=1, RF=2, RR=3
        let orientation = PairOrientation(rawValue: bwaDir)!

        return (orientation, dist)
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

            // Proper-pair bounds for mate rescue and pairing (bwa-mem2 MAPPING_BOUND=3, MAX_STDDEV=4)
            // bwa-mem2 widens to the broader of IQR-based and stddev-based bounds
            let pLow = q25 - 3 * iqr
            let pHigh = q75 + 3 * iqr
            dist.stats[idx].properLow = max(1, min(pLow, Int64(mean - 4.0 * stddev + 0.5)))
            dist.stats[idx].properHigh = max(pHigh, Int64(mean + 4.0 * stddev + 0.5))

            if filtered.count > maxCount {
                maxCount = filtered.count
                primaryIdx = idx
            }
        }

        dist.primaryOrientation = PairOrientation(rawValue: primaryIdx) ?? .fr

        // Post-filter: suppress orientations with less than MIN_DIR_RATIO of the
        // max orientation's count. Matches bwa-mem2 mem_pestat() lines 141-147.
        if maxCount > 0 {
            for d in 0..<4 {
                if !dist.stats[d].failed
                    && Double(dist.stats[d].count) < Double(maxCount) * Self.minDirRatio {
                    dist.stats[d].failed = true
                }
            }
        }

        return dist
    }
}
