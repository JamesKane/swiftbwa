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
- [ ] **End-to-end integration tests** — Unit tests are solid but no full-pipeline tests (FASTQ in → SAM out) validated against bwa-mem2 output.

## Not yet implemented (vs bwa-mem2)

- [ ] **Interleaved FASTQ (`-p`)** — Smart pairing from single interleaved file. Needs read-pair splitting in CLI.
- [ ] **Append FASTQ comment (`-C`)** — bwa-mem2 appends raw FASTQ comment as SAM tags. `ReadSequence.comment` field exists but not wired through output (htslib BAMRecord doesn't support raw tag injection).
- [ ] **Output all alignments (`-a`)** — Emit all above-threshold alignments instead of condensing into XA tag. `flagAll` constant defined but output path not modified.
- [ ] **Manual insert size override (`-I`)** — Bypass estimation with user-provided mean, stddev, max, min.
- [ ] **Fixed batch size (`-K`)** — Process fixed number of input bases per batch for reproducibility across thread counts.
- [ ] **Reference header in XR tag (`-V`)** — Output FASTA header annotation in XR:Z aux tag.
- [ ] **Verbosity levels (`-v`)** — Gate stderr diagnostic output by level (1=error, 2=warning, 3=message, 4+=debug).
- [ ] **Re-seeding (`-y`)** — Additional seeding rounds for long reads with high-occurrence seeds.
- [ ] **Dedup/patch overlapping hits** — `mem_sort_dedup_patch` merges overlapping alignment regions via global re-alignment.
- [ ] **Z-dropoff in SW extension** — `zDrop` is parsed but not applied during Smith-Waterman extension.
- [ ] **Match score scaling validation** — bwa-mem2 only scales when `-A` differs from user defaults; current implementation always scales.
