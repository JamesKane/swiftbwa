# swiftbwa

A Swift 6 reimplementation of [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) with optional Metal GPU acceleration.

swiftbwa builds and reads bwa-mem2-compatible index files and produces equivalent SAM/BAM output, with a concurrency-safe architecture built on Swift's actor model and strict sendability checking.

## Features

- **Native index building**: `swiftbwa index` builds bwa-mem2-compatible indices from FASTA (plain or gzipped) — no need to install bwa-mem2 separately
- **Full bwa-mem2 pipeline**: SMEM seeding, seed chaining, banded Smith-Waterman extension, CIGAR generation
- **Paired-end support**: two-file and interleaved modes, insert size estimation, mate rescue
- **Metal GPU acceleration**: `--gpu` flag offloads Smith-Waterman kernels to Apple Silicon GPU using SIMD wavefront parallelism
- **Streaming I/O**: processes reads in chunks without loading entire FASTQ files into memory
- **Swift 6 concurrency**: actor-based pipeline orchestration, value types everywhere, strict sendability

## Requirements

- macOS 14+
- Swift 6.0+
- [swift-htslib](https://github.com/jkane/swift-htslib) cloned at `../swift-htslib`

## Building

```bash
swift build              # Debug build
swift build -c release   # Release build (enables -Ounchecked on alignment hot paths)
```

## Usage

swiftbwa is CLI-compatible with bwa-mem2. It uses the same index files and accepts the same flags.

### Indexing

```bash
# Build index from FASTA (plain or gzipped)
swiftbwa index ref.fa

# Custom output prefix
swiftbwa index -p /path/to/index ref.fa.gz
```

This produces `.bwt.2bit.64`, `.pac`, `.ann`, and `.amb` files that are byte-identical to bwa-mem2's output. Indices built by either tool are fully interchangeable.

### Alignment

```bash
# Single-end alignment
swiftbwa mem -t 8 ref.fa reads.fq > out.sam

# Paired-end alignment
swiftbwa mem -t 8 ref.fa reads_R1.fq reads_R2.fq > out.sam

# Paired-end, interleaved FASTQ
swiftbwa mem -t 8 -p ref.fa interleaved.fq > out.sam

# BAM output
swiftbwa mem -t 8 -o out.bam ref.fa reads_R1.fq reads_R2.fq

# With Metal GPU acceleration
swiftbwa mem -t 8 --gpu ref.fa reads_R1.fq reads_R2.fq > out.sam
```

### Key alignment options

| Flag | Description |
|------|-------------|
| `-t` | Number of threads (default: 1) |
| `-k` | Minimum seed length (default: 19) |
| `-w` | Band width (default: 100) |
| `-A` | Match score (default: 1) |
| `-B` | Mismatch penalty (default: 4) |
| `-O` | Gap open penalty (default: 6) |
| `-E` | Gap extension penalty (default: 1) |
| `-T` | Minimum alignment score to output (default: 30) |
| `-R` | Read group header line |
| `-M` | Mark shorter split hits as secondary |
| `--gpu` | Enable Metal GPU acceleration |
| `-K` | Batch size in bases for reproducibility |

## Architecture

```
swiftbwa (CLI)
    ├── BWAMem (pipeline orchestration, SAM output)
    │       ├── Alignment (Smith-Waterman, chaining, filtering)
    │       │       └── BWACore
    │       ├── FMIndex (BWT, suffix array, SMEM search)
    │       │       └── BWACore
    │       ├── MetalSW (GPU compute shaders)
    │       │       └── BWACore
    │       └── BWACore (value types: bases, seeds, chains, scoring params)
    └── FMIndexBuilder (FASTA → index construction)
            ├── FMIndex
            └── BWACore
```

- **BWACore** — Zero-dependency value types, all `Sendable`. 2-bit base encoding.
- **FMIndex** — Loads bwa-mem2 index files via mmap. `BWT`, `SuffixArray`, `PackedReference` own raw pointers as `@unchecked Sendable` final classes.
- **Alignment** — Stateless pure functions. Compiled with `-Ounchecked` in release. Includes banded Needleman-Wunsch (raw pointer DP), SIMD16 Farrar striped local SW for mate rescue, seed chaining, and CIGAR generation.
- **MetalSW** — Metal compute shaders for Smith-Waterman. Wavefront kernel uses 32 threads per alignment via SIMD group shuffles for anti-diagonal parallelism.
- **BWAMem** — `BWAMemAligner` actor orchestrates the full pipeline. `SAMOutputBuilder` produces BAM records via swift-htslib.
- **FMIndexBuilder** — Builds bwa-mem2-compatible indices from FASTA. SA-IS suffix array construction, BWT derivation, FM-index checkpoint generation.

## Testing

```bash
swift test                                    # All tests
swift test --filter AlignmentTests            # One suite
swift test --filter "FMIndexTests/testSMEM"   # Single test
```

## License

BSD 3-Clause
