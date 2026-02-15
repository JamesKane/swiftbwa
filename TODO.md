# SwiftBWA — GPU Optimization Opportunities

Analysis of Pham et al., "Accelerating BWA-MEM Read Mapping on GPUs" (ICS '23) against swiftbwa's current Metal implementation.

**Status: All applicable techniques implemented or evaluated. No remaining opportunities from this paper.**

## Implemented

| Paper Technique | swiftbwa Implementation |
|---|---|
| K-mer hash for BWT initialization | 12-mer hash table in `smem_forward.metal` |
| Warp-level SMEM forward extension | 32-thread SIMD groups in `smem_forward.metal` |
| GPU internal reseeding | `internal_reseed.metal`, one thread per task |
| Wavefront anti-diagonal local SW | `local_sw_wavefront.metal`, SIMD shuffles, tile-based (160bp) |
| Bump allocator (CPU side) | `ReadArena` 256KB arena per read |
| Async I/O pipelining | Double-buffered GPU SMEM dispatch |
| GPU global SW + CIGAR | `global_sw.metal` banded NW with traceback |
| Metal buffer pooling | `MetalBufferPool` reusable pool |
| Read reordering for cache coherence | 20-base prefix sort in `SMEMDispatcher` (perf-neutral on M4 Max) |

## Evaluated and skipped

| Paper Technique | Reason |
|---|---|
| GPU-parallel seed chaining | Poor GPU utilization (49% warp efficiency per paper). Seed counts per read too small. CPU chaining with arena already fast. |
| GPU SAM/BAM output generation | <5% of wall time. Would bypass/duplicate swift-htslib for negligible gain. |
| Wavefront early-start for banded SW | GPU banded SW rejected entirely — per-seed tasks are 1-5us, batching overhead exceeds savings. Only local SW (rescue) uses GPU, already wavefront. |
| Read reordering (perf impact) | Implemented but performance-neutral on Apple Silicon. Unified memory eliminates the discrete GPU cache divergence penalty the paper targets (40% on NVIDIA → ~0% on M4 Max). |
