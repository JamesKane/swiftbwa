# SwiftBWA — Remaining Work

## Critical (output correctness)

- [x] **CLI CIGAR fix** — `Sources/swiftbwa/main.swift` has its own hardcoded all-M CIGAR, bypassing `CIGARGenerator`. Needs to use `alignBatch()` or call `CIGARGenerator` directly.
- [x] **MD tag** — Computed in `CIGARGenerator.generateMD()`, written via `aux.updateString(tag: "MD", ...)`.

## Important (real-world usability)

- [x] **Paired-end alignment** — `alignPairedBatch()` wired into CLI. `PairedEndResolver` resolves best pair, `InsertSizeEstimator` infers distribution, TLEN set, proper-pair flagged, MC tag written. Mate rescue implemented (`MateRescue.rescue()` via `LocalSWAligner`) as Phase 2.5 in the paired-end pipeline.
- [x] **CLI parallelization** — `alignBatch()`/`alignPairedBatch()` use TaskGroup with sliding-window concurrency capped at `numThreads`. `alignRead()` is `nonisolated` so child tasks run truly in parallel.
- [x] **Supplementary/chimeric alignments** — Two-pass output classification (primary/supplementary/secondary) in `BWAMemAligner.emitSingleEndAlignments()`. SA and XA tag generation. Hard-clip conversion for supplementary records. MAPQ capping. Paired-end path emits supplementary alignments with proper mate info.

## Nice-to-have / advanced

- [x] **ALT-aware alignment** — `.alt` file loading in `FMIndexLoader`, `isAlt` propagated through chain→region pipeline, two-phase `markSecondaryALT()` prevents ALT hits from suppressing primaries, `pa:f` tag output, ALT-aware XA limits (200), MAPQ cap skip for ALT supplementary hits.
- [x] **Read group support** — `-R` flag accepts `@RG` header line, inserted verbatim via `addLines()`. `readGroupID` extracted and written as `RG` aux tag on all output records.
- [x] **Clipping penalties** — `penClip5`/`penClip3` used in `ExtensionAligner.extend()` clip-vs-extend decision (bwa-mem2 logic). Creates split-read candidates when local alignment beats end-to-end minus clip penalty.
- [x] **CLI flags** — Full bwa-mem2-compatible flag set: `-A` (match score + scaling), `-B`, `-O`/`-E` (separate ins/del), `-L` (clip penalties), `-U` (unpaired penalty), `-r` (seed split ratio), `-c` (max occurrences), `-D` (chain drop ratio), `-W` (min chain weight), `-m` (max mate rescue), `-h` (XA limits), `-M` (mark splits as secondary), `-Y` (soft-clip supplementary), `-5` (primary5 reorder), `-q` (keep supp MAPQ), `-S` (skip rescue), `-P` (skip pairing), `-j` (ignore ALT), `-H` (custom header lines).
- [x] **End-to-end integration tests** — 10 E2E tests (7 single-end, 3 paired-end) validate full pipeline (FASTQ in → SAM out) against bwa-mem2 gold-standard output on lambda phage. Uncovered and fixed two bugs: compressed SA resolution (BWT walk-back) and reverse-strand FLAG logic.

## Not yet implemented (vs bwa-mem2)

- [x] **Interleaved FASTQ (`-p`)** — `-p` flag deinterleaves single FASTQ into read1/read2 by index parity, feeds into existing `alignPairedBatch()`.
- [x] **Append FASTQ comment (`-C`)** — FASTQ comment preserved via `readFASTQWithComments()` line-by-line reader (handles gzip), written as `CO:Z` aux tag via `SAMOutputBuilder`.
- [x] **Output all alignments (`-a`)** — Emit all above-threshold alignments as full SAM records with 0x100 flag instead of condensing into XA tag. XA tag skipped when `-a` is set.
- [x] **Manual insert size override (`-I`)** — `InsertSizeOverride` struct with bwa-mem2 defaults for omitted max/min. `InsertSizeEstimator.buildManualDistribution()` constructs FR-orientation stats.
- [x] **Fixed batch size (`-K`)** — `alignSEInChunks`/`alignPairedInChunks` process reads in chunks of K input bases. Insert size estimated per chunk.
- [x] **Reference header in XR tag (`-V`)** — `outputRefHeader` option writes `XR:Z:<anno>` aux tag from `ReferenceAnnotation.anno`.
- [x] **Verbosity levels (`-v`)** — `options.verbosity` gates all `fputs` calls: 1=error, 2=warning, 3=info (default), 4+=debug.
- [x] **Re-seeding (`-y`)** — Phase 1.5 in `alignRead()`: re-runs SMEM finding with `minIntv=maxOccurrences` for high-occurrence seeds, merges/deduplicates results.
- [x] **Dedup/patch overlapping hits** — `mem_sort_dedup_patch` merges overlapping alignment regions via global re-alignment. Implemented as `RegionDedup.sortDedupPatch()` in Phase 4.5 of the alignment pipeline.
- [x] **Z-dropoff in SW extension** — SIMD implementations (BandedSW8/BandedSW16) rewritten with proper Farrar striped algorithm: diagonal shift, lazy-F correction, h0 initialization, z-dropoff, globalScore/globalTargetEnd/maxOff tracking. ExtensionAligner now uses tiered SIMD (8-bit → 16-bit fallback).
- [x] **Match score scaling validation** — Penalty scaling now guarded by `if a != 1`, matching bwa-mem2 behavior.
