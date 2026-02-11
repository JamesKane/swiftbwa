# SwiftBWA — Remaining Work

## Critical (output correctness)

- [x] **CLI CIGAR fix** — `Sources/swiftbwa/main.swift` has its own hardcoded all-M CIGAR, bypassing `CIGARGenerator`. Needs to use `alignBatch()` or call `CIGARGenerator` directly.
- [ ] **MD tag** — Not implemented. Required by GATK, samtools calmd, most variant callers. Walk CIGAR + reference bases to produce MD string.

## Important (real-world usability)

- [ ] **Paired-end alignment** — `PairedEndResolver` and `InsertSizeEstimator` exist but are unwired. Missing: mate rescue, insert size inference, TLEN, proper pair orientation, MC tag.
- [ ] **CLI parallelization** — Actor + TaskGroup infrastructure exists in `alignBatch()` but CLI processes reads serially. Wire up batching and honor `-t` thread flag.
- [ ] **Supplementary/chimeric alignments** — No SA or XA tag generation. All non-primary alignments are secondary, never supplementary. No split-read detection.

## Nice-to-have / advanced

- [ ] **ALT-aware alignment** — `MemAlnReg.isAlt` field exists but always `false`. No `.alt` file loading or ALT-contig logic.
- [ ] **Read group support** — No `@RG` header line or RG tags on records.
- [ ] **Clipping penalties** — `penClip5`/`penClip3` in `ScoringParameters` are unused in scoring decisions.
- [ ] **End-to-end integration tests** — Unit tests are solid but no full-pipeline tests (FASTQ in → SAM out) validated against bwa-mem2 output.
