#include <metal_stdlib>
using namespace metal;

/// GPU forward extension for SMEM finding.
/// One SIMD group (32 threads) per read. Each thread extends rightward from
/// its assigned query position using FM-index backward search on BWT checkpoints.
///
/// K-mer hash table optimization: pre-computed BWT intervals for all 4^12 possible
/// 12-mers. Each thread does an O(1) table lookup to skip the first 11 BWT extensions.
/// Positions where the 12-mer has no match (interval size 0) are skipped entirely
/// since they can't produce a valid SMEM (minSeedLen is typically 19 > 12).
///
/// Results per position: (rightEnd, k, s) where rightEnd is the exclusive right endpoint,
/// k is the forward SA interval start, and s is the interval size.
/// The caller identifies left-maximal positions on CPU by comparing adjacent rightEnds.

#define SIMD_W 32
#define KMER_K 12

/// Checkpoint occurrence structure — must match Swift's CheckpointOCC layout exactly.
/// 4 × int64 counts + 4 × uint64 bitstrings = 64 bytes.
struct CheckpointOCC {
    int64_t counts[4];      // cumulative A, C, G, T counts at this checkpoint
    uint64_t bitstrings[4]; // one-hot encoded BWT characters for 64 positions
};

/// Pre-computed BWT interval for a 12-mer. 24 bytes per entry.
struct KmerEntry {
    int64_t k;  // forward SA interval start
    int64_t l;  // reverse SA interval start
    int64_t s;  // interval size
};

/// Compute one-hot mask for rank query at position y within a 64-char checkpoint.
inline uint64_t one_hot_mask(int y) {
    if (y <= 0) return 0;
    if (y >= 64) return 0xFFFFFFFFFFFFFFFF;
    return ~(uint64_t(0xFFFFFFFFFFFFFFFF) >> uint(y));
}

/// Get occurrence count of base c (0-3) at position pos in the BWT.
inline int64_t occ(device const CheckpointOCC* cps, int64_t pos, int c) {
    int occID = int(pos >> 6);  // pos / 64
    int y = int(pos & 63);      // pos % 64

    int64_t cp_count = cps[occID].counts[c];
    uint64_t bits = cps[occID].bitstrings[c];
    uint64_t mask = one_hot_mask(y);
    return cp_count + int64_t(popcount(bits & mask));
}

/// Perform backward extension on a bidirectional SA interval.
/// Computes all 4 bases' k/s values and the l ladder.
inline void backward_ext(
    device const CheckpointOCC* cps,
    device const int64_t* bwt_counts,
    int64_t sentinel,
    int64_t k_in, int64_t l_in, int64_t s_in,
    int a,
    thread int64_t& k_out, thread int64_t& l_out, thread int64_t& s_out
) {
    int64_t sp = k_in;
    int64_t ep = k_in + s_in;

    int64_t kAll[4], sAll[4];
    for (int c = 0; c < 4; c++) {
        int64_t occ_sp = occ(cps, sp, c);
        int64_t occ_ep = occ(cps, ep, c);
        kAll[c] = bwt_counts[c] + occ_sp;
        sAll[c] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = (sp <= sentinel && ep > sentinel) ? 1 : 0;
    int64_t lAll[4];
    lAll[3] = l_in + sentinel_offset;
    lAll[2] = lAll[3] + sAll[3];
    lAll[1] = lAll[2] + sAll[2];
    lAll[0] = lAll[1] + sAll[1];

    k_out = kAll[a];
    l_out = lAll[a];
    s_out = sAll[a];
}

/// Forward extension step: extend the match rightward by one base.
inline bool forward_ext(
    device const CheckpointOCC* cps,
    device const int64_t* bwt_counts,
    int64_t sentinel,
    int base,
    thread int64_t& k, thread int64_t& l, thread int64_t& s
) {
    int64_t swapped_k = l;
    int64_t swapped_l = k;
    int64_t swapped_s = s;

    int64_t new_k, new_l, new_s;
    backward_ext(cps, bwt_counts, sentinel,
                 swapped_k, swapped_l, swapped_s,
                 3 - base,
                 new_k, new_l, new_s);

    k = new_l;
    l = new_k;
    s = new_s;

    return s > 0;
}

kernel void smem_forward(
    device const uchar*          queries        [[buffer(0)]],
    device const uint*           queryOffsets   [[buffer(1)]],
    device const ushort*         queryLengths   [[buffer(2)]],
    device const CheckpointOCC*  checkpoints    [[buffer(3)]],
    device const int64_t*        bwtCounts      [[buffer(4)]],  // 5 values
    device const int64_t*        sentinel       [[buffer(5)]],  // 1 value
    device       int*            rightEnds      [[buffer(6)]],  // per-position
    device       int64_t*        kValues        [[buffer(7)]],  // per-position
    device       int64_t*        sValues        [[buffer(8)]],  // per-position
    device const uint*           resultOffsets  [[buffer(9)]],  // per-read
    device const uint*           numTasks       [[buffer(10)]],
    device const KmerEntry*      kmerTable      [[buffer(11)]],  // 4^12 entries
    uint lane    [[thread_index_in_simdgroup]],
    uint sg_idx  [[simdgroup_index_in_threadgroup]],
    uint tg_idx  [[threadgroup_position_in_grid]],
    uint sg_per_tg [[simdgroups_per_threadgroup]]
) {
    uint task = tg_idx * sg_per_tg + sg_idx;
    if (task >= numTasks[0]) return;

    uint qlen = queryLengths[task];
    device const uchar* Q = queries + queryOffsets[task];
    uint resOff = resultOffsets[task];
    int64_t sent = sentinel[0];

    // Process positions in rounds of 32 (SIMD width)
    uint nRounds = (qlen + SIMD_W - 1) / SIMD_W;

    for (uint round = 0; round < nRounds; round++) {
        uint qi = round * SIMD_W + lane;

        if (qi >= qlen) {
            continue;
        }

        int64_t cur_k, cur_l, cur_s;
        int rEnd;
        uint startJ;
        bool initialized = false;

        // Try K-mer hash lookup: skip 11 BWT extensions with O(1) table read
        if (qi + KMER_K <= qlen) {
            int hashValue = 0;
            bool validKmer = true;
            for (int i = 0; i < KMER_K; i++) {
                int b = int(Q[qi + i]);
                if (b >= 4) { validKmer = false; break; }
                hashValue = hashValue * 4 + b;
            }

            if (validKmer) {
                KmerEntry entry = kmerTable[hashValue];
                if (entry.s > 0) {
                    // Fast path: use pre-computed interval, extend from qi+12
                    cur_k = entry.k;
                    cur_l = entry.l;
                    cur_s = entry.s;
                    rEnd = int(qi + KMER_K);
                    startJ = qi + KMER_K;
                    initialized = true;
                } else {
                    // 12-mer has no match → no SMEM >= 12bp starting here
                    rightEnds[resOff + qi] = int(qi);
                    kValues[resOff + qi] = 0;
                    sValues[resOff + qi] = 0;
                    continue;
                }
            }
            // If !validKmer (N in first 12 bases), fall through to single-base init
        }

        if (!initialized) {
            // Fallback: near end of read or N in K-mer window
            int base = int(Q[qi]);
            if (base >= 4) {
                rightEnds[resOff + qi] = int(qi);
                kValues[resOff + qi] = 0;
                sValues[resOff + qi] = 0;
                continue;
            }

            cur_k = bwtCounts[base];
            cur_l = bwtCounts[3 - base];
            cur_s = bwtCounts[base + 1] - bwtCounts[base];

            if (cur_s <= 0) {
                rightEnds[resOff + qi] = int(qi);
                kValues[resOff + qi] = 0;
                sValues[resOff + qi] = 0;
                continue;
            }

            rEnd = int(qi + 1);
            startJ = qi + 1;
        }

        // Extend rightward from startJ to end of read
        for (uint j = startJ; j < qlen; j++) {
            int b = int(Q[j]);
            if (b >= 4) break;

            int64_t new_k = cur_k, new_l = cur_l, new_s = cur_s;
            bool ok = forward_ext(checkpoints, bwtCounts, sent, b,
                                  new_k, new_l, new_s);

            if (!ok || new_s <= 0) break;

            cur_k = new_k;
            cur_l = new_l;
            cur_s = new_s;
            rEnd = int(j + 1);
        }

        rightEnds[resOff + qi] = rEnd;
        kValues[resOff + qi] = cur_k;
        sValues[resOff + qi] = cur_s;
    }
}
