// swift-tools-version: 6.0
import PackageDescription

let package = Package(
    name: "swiftbwa",
    platforms: [.macOS(.v14)],
    products: [
        .executable(name: "swiftbwa", targets: ["swiftbwa"]),
        .library(name: "BWAMem", targets: ["BWAMem"]),
        .library(name: "Alignment", targets: ["Alignment"]),
        .library(name: "FMIndex", targets: ["FMIndex"]),
        .library(name: "BWACore", targets: ["BWACore"]),
    ],
    dependencies: [
        .package(path: "../swift-htslib"),
        .package(url: "https://github.com/apple/swift-argument-parser.git", from: "1.3.0"),
    ],
    targets: [
        // Core types — zero dependencies
        .target(
            name: "BWACore",
            swiftSettings: [.enableExperimentalFeature("StrictConcurrency")]
        ),

        // FM-Index loading and search
        .target(
            name: "FMIndex",
            dependencies: ["BWACore"],
            swiftSettings: [.enableExperimentalFeature("StrictConcurrency")]
        ),

        // Alignment: Smith-Waterman, chaining, filtering
        .target(
            name: "Alignment",
            dependencies: ["BWACore"],
            swiftSettings: [.enableExperimentalFeature("StrictConcurrency")]
        ),

        // BWAMem: pipeline orchestration + SAM output
        .target(
            name: "BWAMem",
            dependencies: [
                "BWACore",
                "FMIndex",
                "Alignment",
                .product(name: "Htslib", package: "swift-htslib"),
            ],
            swiftSettings: [.enableExperimentalFeature("StrictConcurrency")]
        ),

        // CLI executable
        .executableTarget(
            name: "swiftbwa",
            dependencies: [
                "BWAMem",
                .product(name: "ArgumentParser", package: "swift-argument-parser"),
            ],
            swiftSettings: [.enableExperimentalFeature("StrictConcurrency")]
        ),

        // Tests
        .testTarget(name: "BWACoreTests", dependencies: ["BWACore"]),
        .testTarget(name: "FMIndexTests", dependencies: ["FMIndex"]),
        .testTarget(name: "AlignmentTests", dependencies: ["Alignment"]),
        .testTarget(name: "BWAMemTests", dependencies: ["BWAMem", "FMIndex", "Alignment"]),
    ]
)
