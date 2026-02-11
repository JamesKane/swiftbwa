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
- [ ] **End-to-end integration tests** — Unit tests are solid but no full-pipeline tests (FASTQ in → SAM out) validated against bwa-mem2 output.
