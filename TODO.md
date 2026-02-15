# SwiftBWA — ML/NN Acceleration Roadmap

## Literature Review Summary

Analysis of neural network and machine learning techniques applicable to BWA-MEM-style alignment pipelines, evaluated for Apple Silicon (M4 Max) integration.

### Sources

| Paper | Year | Venue | Targets |
|-------|------|-------|---------|
| BWA-MEME (Jung & Han) | 2022 | Bioinformatics | Learned index for SA lookup |
| Sapling (Kirsche et al.) | 2021 | Bioinformatics | Piecewise linear SA predictor |
| LISA (Ho et al.) | 2021 | bioRxiv | IP-BWT + RMI for SMEM search |
| Cline et al. | 2020 | PeerJ | ML-based MAPQ recalibration |
| GateKeeper-GPU (Yaman et al.) | 2024 | IEEE TPDS | Pre-alignment Hamming filter |
| AGNES (Arafat et al.) | 2025 | arXiv | GNN seed chaining (long reads) |
| ESA / dna2vec (Holur et al.) | 2025 | Bioinformatics | Transformer alignment |
| CARE 2.0 (Kallenborn et al.) | 2022 | BMC Bioinformatics | Read error correction |
| E-score (Ashrafzadeh et al.) | 2024 | Briefings in Bioinformatics | Embedding-based scoring (protein) |
| Differentiable SW (Koo et al.) | 2023 | Bioinformatics | Differentiable alignment (training) |

### Evaluation Criteria

- **Latency budget**: Per-read operations must complete in microseconds. Per-batch operations can tolerate milliseconds.
- **Accuracy**: Must not degrade alignment quality vs. bwa-mem2 baseline.
- **Apple Silicon fit**: Prefer Accelerate `vDSP`/BNNS or pure Swift for small models. CoreML/ANE dispatch latency (~1ms) rules out per-read inference.
- **Bottleneck relevance**: Current profile shows 43.5% memory management, ~11% algorithmic work (6.1% banded SW, 4.6% BWT/SMEM). MAPQ accuracy gaps remain (91.4% vs 98.1% proper pair rate).

### Rejected Approaches

| Approach | Reason |
|----------|--------|
| Transformer alignment (ESA) | 100x slower than conventional SW |
| GNN seed chaining (AGNES) | Long-read only, 1ms/read overhead, synthetic validation only |
| Differentiable SW | Training tool, inherently slower than standard SW |
| CoreML/ANE per-read inference | ~1ms dispatch latency, 1000x too slow |
| DNA foundation model embeddings | Wrong abstraction level for alignment inner loops |
| Inline error correction (CARE 2.0) | Requires whole-dataset pre-processing, not inline |
| LISA learned index | 33-92 GB memory overhead, designed for Intel servers |

---

## Experiment Roadmap

### Experiment 1: Learned MAPQ Model

**Goal**: Replace hand-tuned MAPQ formula with a small decision tree ensemble to close accuracy gaps.

**Why**: Known remaining gaps — 91.4% vs 98.1% proper pair rate, 98.2% vs 99.4% mapping rate, ~20K overconfident MAPQ at 1M reads. A learned model trained on bwa-mem2 outputs should capture the complex interactions between alignment features that the heuristic formula approximates poorly.

**Approach**:
1. Generate training data: align 10M+ reads with both swiftbwa and bwa-mem2, extract per-read feature vectors + bwa-mem2 MAPQ as labels
2. Features (all already computed in pipeline): alignment score, sub-optimal score (`sub`), number of chains, seed count, edit distance, CIGAR complexity, paired-end insert size deviation, read length, number of Ns
3. Train gradient-boosted decision tree (XGBoost/LightGBM) on MAPQ regression or classification into MAPQ bins
4. Export as flat threshold arrays — nested if-else in pure Swift
5. Evaluate: compare proper pair rate, mapping rate, MAPQ distribution, SNP calling concordance

**Implementation**: Pure Swift, no framework dependencies. <100ns per read.

**Validation**:
- MAPQ distribution comparison vs bwa-mem2 on held-out reads
- Variant calling concordance (GATK HaplotypeCaller) on HG002 truth set
- No regression in alignment positions or CIGAR strings (MAPQ-only change)

**Risk**: Low — MAPQ is a post-hoc quality annotation. Wrong MAPQ doesn't change where reads map, only downstream filtering thresholds. Can A/B test against current formula.

---

### Experiment 2: P-RMI Learned Index for Suffix Array Lookup

**Goal**: Replace binary search over suffix array with learned model inference + bounded search, following BWA-MEME.

**Why**: BWA-MEME reports 3.45x seeding throughput and 1.42x end-to-end speedup with identical output. Seeding is currently 4.6% of CPU time (~13.5s CPU on 4M PE reads), so absolute gain is bounded but real.

**Approach**:
1. Port BWA-MEME's P-RMI training (Rust → Swift or Python script)
   - 3-layer recursive model: layer 0 = 1 linear model, layer 1 = N linear models, layer 2 = M linear models
   - Each model: `predicted_pos = slope * key + intercept` (2 FP64 values)
   - Training: sort suffix array entries, fit piecewise linear models with bounded error
2. Index build: train P-RMI on suffix array, store alongside existing `.bwt.2bit.64` index
3. Inference: encode query k-mer → FP64 key → P-RMI predicts SA position → bounded binary search within error bound
4. Integrate into `SuffixArray.lookup()` path — transparent to `SMEMFinder`

**Implementation**: Accelerate `vDSP_vmaD` for batch multiply-accumulate, or plain Swift FP64 arithmetic. Model weights are a few MB (fits in L2 cache).

**Validation**:
- Bit-identical SA lookup results (the model only narrows the search, doesn't change the answer)
- Benchmark seeding phase isolation and end-to-end on 4M PE reads
- Memory overhead measurement (model size vs. existing index)

**Risk**: Low — produces identical results by construction. The bounded search guarantees correctness. Main risk is that the speedup may be smaller on Apple Silicon (unified memory, different cache hierarchy than Intel servers where BWA-MEME was benchmarked).

---

### Experiment 3: Piecewise Linear SA Predictor (Sapling)

**Goal**: Simpler alternative to P-RMI — a single piecewise linear function predicting SA positions.

**Why**: <1% memory overhead (vs P-RMI's multi-MB model), ~2x seeding speedup, trivial to implement. Good low-risk experiment to establish whether learned SA prediction helps at all on Apple Silicon before investing in the full P-RMI.

**Approach**:
1. Build piecewise linear model: divide suffix array into N segments, fit slope + intercept per segment
2. At query time: hash query → select segment → predict position → binary search within error bound
3. Integrate into `SuffixArray` as an optional acceleration path

**Implementation**: Pure Swift — array of `(slope: Double, intercept: Double)` pairs, one comparison + one multiply-add per lookup. No framework dependencies.

**Validation**: Same as Experiment 2 — bit-identical results, benchmark seeding phase.

**Risk**: Minimal. If speedup is negligible on Apple Silicon, the code is small enough to discard.

---

### Experiment 4: GPU Pre-Alignment Hamming Filter

**Goal**: Discard clearly-bad seed extension candidates before expensive banded SW, using a fast GPU Hamming distance kernel.

**Why**: GateKeeper-GPU reports 2.9x alignment acceleration by filtering 52x fewer false accepts. If many candidate extension regions are being sent to banded SW that ultimately fail the score threshold, a cheap pre-filter saves expensive DP work.

**Approach**:
1. Measure current reject rate: count how many `ExtensionAligner` calls produce alignments below `T` (score threshold 30)
2. If reject rate is high (>20%), implement a Metal kernel computing Hamming distance between query and candidate reference region
3. Filter candidates exceeding a Hamming distance threshold before dispatching banded SW
4. Integrate into the extension phase, between chain selection and SW dispatch

**Implementation**: Custom Metal compute kernel (one thread per candidate), or `MPSGraph` Hamming distance operation.

**Prerequisite**: Profiling data showing significant time spent on ultimately-rejected extensions. If reject rate is low, this experiment has no payoff.

**Validation**:
- No change in alignment output (filter threshold must be conservative enough to never discard true positives)
- Benchmark extension phase isolation and end-to-end
- False negative analysis: ensure no alignments above `T` are filtered

**Risk**: Medium — requires careful threshold tuning. Too aggressive = missed alignments. Too conservative = no speedup. The reject rate may already be low given chain filtering and score thresholds.

---

## Experiment Priority and Dependencies

```
Experiment 1 (Learned MAPQ)          ← highest value, addresses known accuracy gap
    │
    ├── independent of all others
    │
Experiment 3 (Sapling linear SA)     ← lowest effort, establishes baseline for learned SA
    │
    └── if promising → Experiment 2 (P-RMI full learned index)

Experiment 4 (Hamming pre-filter)    ← conditional on reject rate profiling
```

Experiment 1 should run first — it targets accuracy (the known gap), not performance. Experiments 2-3 target seeding speed (4.6% of CPU time). Experiment 4 targets extension speed but may have no payoff if reject rates are already low.

## Completed: GPU Optimization (Pham et al. ICS '23)

All applicable GPU techniques from Pham et al. have been implemented or evaluated:

| Technique | Status |
|-----------|--------|
| K-mer hash for BWT initialization | Implemented (`smem_forward.metal`) |
| Warp-level SMEM forward extension | Implemented (32-thread SIMD groups) |
| GPU internal reseeding | Implemented (`internal_reseed.metal`) |
| Wavefront anti-diagonal local SW | Implemented (`local_sw_wavefront.metal`) |
| Bump allocator (CPU side) | Implemented (`ReadArena`) |
| Async I/O pipelining | Implemented (double-buffered dispatch) |
| GPU global SW + CIGAR | Implemented (`global_sw.metal`) |
| Metal buffer pooling | Implemented (`MetalBufferPool`) |
| Read reordering for cache coherence | Implemented (perf-neutral on M4 Max) |
| GPU-parallel seed chaining | Skipped (poor GPU utilization, 49% warp efficiency) |
| GPU SAM/BAM output | Skipped (<5% wall time) |
| Wavefront early-start banded SW | Skipped (GPU banded SW rejected, per-seed tasks too small) |
