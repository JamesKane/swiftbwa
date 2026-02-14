# GlobalAligner Performance Optimization — Findings Summary

## Problem

Instruments profiling of 100k SE reads showed **94% of alignment time** (~4047 of 4223 samples) in `GlobalAligner.align()` — the scalar banded Needleman-Wunsch DP for CIGAR generation. This caused a ~4.35x SE slowdown vs bwa-mem2.

## Changes Made

### 1. Raw pointer allocation (GlobalAligner.swift)

Replaced Swift Arrays with `UnsafeMutablePointer.allocate()` for all DP arrays (`qp`, `h`, `e`, `z`). Key findings:

- `UnsafeMutableBufferPointer.subscript` still calls `_checkIndex` → `_precondition` in `-O` builds — NOT unchecked
- `UnsafeMutablePointer.subscript` (raw pointer) IS truly unchecked
- Impact: ~25% improvement in per-call DP cost (16s → 12s at 100k serial)

### 2. Wrapping arithmetic (GlobalAligner.swift)

Changed all `Int32` arithmetic in the hot inner loop to wrapping operators (`&+`, `&-`, `&*`). Standard operators trap on overflow even in `-O` builds. This is required to match C's behavior.

### 3. Branchless ternary restructuring (GlobalAligner.swift)

Restructured `if/else` blocks to ternary expressions + `max()` to encourage ARM64 `CSEL` (conditional select) code generation, matching bwa-mem2's `ksw_global2()` style. Minimal additional impact (~0.1s).

### 4. Pre-computed scoring matrix (all 3 files)

- `CIGARGenerator.generate()` and `BWAMemAligner.generateCIGAR()` accept optional `scoringMatrix: [Int8]?`
- Matrix built once per batch, passed through call chain
- Eliminates ~100k `malloc`/`free` cycles per batch

### 5. Parallel SE CIGAR generation (BWAMemAligner.swift) — **BIGGEST WIN**

SE path was running CIGAR generation serially in `emitSingleEndAlignments()`. Restructured `alignBatch()` to combine alignment + CIGAR generation in one parallel task:

```swift
group.addTask { [self] in
    let regions = self.alignRead(read, readId: UInt64(idx))
    let mat = self.options.scoring.scoringMatrix()
    let cigars = regions.map { self.generateCIGAR(read: read, region: $0, scoringMatrix: mat) }
    return (idx, regions, cigars)
}
```

Impact: SE 100k dropped from ~12s → ~4.4s (serial CIGAR was the actual bottleneck, not per-call DP speed).

### 6. Merged PE rescue + CIGAR phases (BWAMemAligner.swift)

Combined the separate rescue and CIGAR pre-computation task groups into one parallel phase. This eliminates a synchronization barrier between the two phases.

## Results

| Tier | Mode | Before | After | bwa-mem2 | Ratio | Status |
|------|------|--------|-------|----------|-------|--------|
| 100k | SE | ~16s | 4.4s | 3.6s | 1.22x | PASS |
| 100k | PE | ~18s | 17.9s | 9.3s | 1.92x | PASS |
| 1M | SE | — | 28.6s | 15.2s | 1.88x | PASS |
| 1M | PE | — | 170.2s | 71.1s | 2.39x | **FAIL** |

All 160 tests pass. Threshold is 2.0x bwa-mem2.

## Root Cause of Remaining PE 1M Gap

The per-call `GlobalAligner.align()` is still **~7x slower** than bwa-mem2's C `ksw_global2()`. Swift-level optimizations (raw pointers, wrapping arithmetic, branchless ternaries) recovered ~25% but cannot close the compiler quality gap. At 1M PE reads, the sheer volume of CIGAR calls (2x reads × multiple regions × rescue) makes this the dominant cost.

bwa-mem2's `ksw_global2()` (ksw.cpp:561-671) advantages:
- C compiler produces tighter inner loop code (fewer register spills, better instruction scheduling)
- Interleaved `eh_t` struct (H and E in one cache line) vs separate `h`/`e` arrays
- `int8_t` query profile (1 byte) vs `Int32` (4 bytes) — 4x better cache utilization
- No function call overhead for pointer arithmetic

## Recommended Next Step

Write the DP inner loop in C as a new SPM target and call it from Swift:

```
Sources/AlignmentCore/          (C target)
  include/alignment_core.h      (public header)
  global_align_dp.c             (DP fill + query profile + backpointer matrix)
```

The C function would handle:
- Query profile building (int8_t)
- DP fill with interleaved eh_t struct
- Backpointer matrix population

Traceback and CIGAR building would stay in Swift (they're not hot — only ~30 samples in profiling).

This should close the ~7x per-call gap and bring PE 1M under the 2x threshold.

---

## Internal Seeding for Long SMEMs (Feb 2026)

### Problem

Position comparison between swiftbwa and bwa-mem2 (1M PE reads, chm13v2.0) revealed
55,798 reads where swiftbwa assigned MAPQ >= 50 but bwa-mem2 assigned MAPQ <= 5 at
different positions. 73% of these had identical alignment scores — both tools found
equally good alignments, but bwa-mem2 found *multiple* while swiftbwa found only one.

**Root cause**: bwa-mem2 performs "internal seeding" (`mem_collect_intv`, bwamem.cpp:695-753).
For long SMEMs (>= minSeedLen x 1.5) with low occurrence (<= 10), it re-searches from the
SMEM midpoint with `min_intv = occurrence + 1` to discover shorter seeds at alternative
genomic loci. swiftbwa had the parameters (`seedSplitRatio`, `splitWidth`) but never used
them.

### Fix

Added `internalReseed()` method to `BWAMemAligner` that splits long low-occurrence SMEMs
at their midpoints. Called from both `alignReadPhase1to3` (CPU path) and
`alignReadPhase1_5to3` (GPU SMEM path).

Also fixed `seedLen0` bug in `ExtensionAligner` — was set to `seeds.count` (number of
seeds in chain) instead of `seed.len` (individual seed length).

### Results (1M PE reads, chm13v2.0, M4 Max, 12 threads)

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| Wall time | 29.6s | 34.9s | +18% |
| CPU time | 425.5s | 499.7s | +17% |
| Overconfident MAPQ (swift>=50, bwa<=5, diff pos) | 55,798 | 38,225 | **-31%** |
| Mean MAPQ for mismatched reads | 32.7 | 25.5 | **-22%** |
| Both-MAPQ-0 agreement (mismatches) | ~few | 54,661 | **much better** |
| Supplementaries for mismatches | ~50K | 35,862 | **-28%** |
| Position concordance | 92.5% | 92.06% | -0.4% |
| High-confidence mismatches (both MAPQ>=30) | ~8.5K | 8,598 | ~same |

The slight position concordance decrease is expected — internal seeding correctly identifies
more multi-mapping reads (swiftbwa MAPQ=0 rose from ~5K to 66K for mismatches), causing
position changes. Reads that move now carry low MAPQ (uncertainty acknowledged) instead of
MAPQ=60 (false confidence).

### What didn't work: Cross-chain seed coverage

bwa-mem2 checks each seed against ALL previously extended regions across all chains before
extending, using query+reference containment and diagonal proximity. We attempted to
replicate this by passing accumulated regions to `ExtensionAligner.extend()`.

This made things **significantly worse** (108K overconfident vs 38K). The fundamental issue:
bwa-mem2 processes all seeds in a single flat loop where sub/subN updates propagate to
covering regions immediately. Our per-chain architecture processes chains sequentially —
the cross-chain coverage suppressed extensions without properly propagating multi-mapping
evidence through sub scores.

### Remaining gaps vs bwa-mem2

- 38K overconfident MAPQ cases remain (was 56K)
- Supplementary count still higher than bwa-mem2 (36K vs 1.5K for mismatches)
- Proper pair rate: ~92.4% vs bwa-mem2's 98.1%
- 18% wall time overhead from additional BWT queries

### Benchmark history

| Configuration | Wall | CPU | vs bwa-mem2 |
|---------------|------|-----|-------------|
| bwa-mem2 2.2.1 | 77.7s | 586.8s | 1.00x |
| swiftbwa CPU (before) | 29.6s | 425.5s | 2.62x faster |
| swiftbwa CPU (internal seeding) | 34.9s | 499.7s | 2.23x faster |
| swiftbwa GPU (before) | 17.4s | 114.0s | 4.47x faster |
