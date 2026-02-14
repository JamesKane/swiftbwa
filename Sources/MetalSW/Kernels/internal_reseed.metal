#include <metal_stdlib>
using namespace metal;

/// GPU kernel for internal reseeding: finds shorter SMEMs at midpoints of
/// long, low-occurrence SMEMs. Faithful reimplementation of
/// SMEMFinder.findSMEMsAtPosition() — one thread per reseed task.

// --- Shared BWT helpers (may already be defined via source concatenation) ---
#ifndef BWT_HELPERS_DEFINED
#define BWT_HELPERS_DEFINED

struct CheckpointOCC {
    int64_t counts[4];
    uint64_t bitstrings[4];
};

inline uint64_t one_hot_mask(int y) {
    if (y <= 0) return 0;
    if (y >= 64) return 0xFFFFFFFFFFFFFFFF;
    return ~(uint64_t(0xFFFFFFFFFFFFFFFF) >> uint(y));
}

inline int64_t occ(device const CheckpointOCC* cps, int64_t pos, int c) {
    int occID = int(pos >> 6);
    int y = int(pos & 63);
    int64_t cp_count = cps[occID].counts[c];
    uint64_t bits = cps[occID].bitstrings[c];
    uint64_t mask = one_hot_mask(y);
    return cp_count + int64_t(popcount(bits & mask));
}

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

#endif // BWT_HELPERS_DEFINED

// --- Internal reseed kernel ---

#define MAX_OUT  16
#define MAX_PREV 96

/// Bidirectional SA interval for thread-local tracking during search.
struct BidiInterval {
    int64_t k;   // Forward SA interval start
    int64_t l;   // Reverse SA interval position
    int64_t s;   // Interval size
    int m;       // Query begin (inclusive)
    int n;       // Query end (inclusive, last matched position)
};

/// One thread per reseed task. Implements findSMEMsAtPosition from a given
/// startPos with minIntv > 1 to find shorter, more specific seeds.
kernel void internal_reseed(
    device const uchar*          queries        [[buffer(0)]],   // packed read bases
    device const uint*           queryOffsets   [[buffer(1)]],   // per-read offset
    device const ushort*         queryLengths   [[buffer(2)]],   // per-read length
    device const uint*           taskReadIdx    [[buffer(3)]],   // which read per task
    device const ushort*         taskStartPos   [[buffer(4)]],   // startPos per task
    device const int64_t*        taskMinIntv    [[buffer(5)]],   // minIntv per task
    device const CheckpointOCC*  checkpoints    [[buffer(6)]],   // BWT checkpoints
    device const int64_t*        bwtCounts      [[buffer(7)]],   // cumulative counts
    device const int64_t*        sentinel       [[buffer(8)]],   // sentinel position
    device       int64_t*        outK           [[buffer(9)]],   // output k (MAX_OUT/task)
    device       int64_t*        outL           [[buffer(10)]],  // output l
    device       short*          outQBegin      [[buffer(11)]],  // output queryBegin
    device       short*          outQEnd        [[buffer(12)]],  // output queryEnd
    device       uchar*          outCount       [[buffer(13)]],  // SMEMs produced/task
    device const uint*           numTasks       [[buffer(14)]],  // total task count
    device const short*          minSeedLenBuf  [[buffer(15)]],  // min seed length (1 val)
    uint tid [[thread_position_in_grid]]
) {
    if (tid >= numTasks[0]) return;

    uint readIdx = taskReadIdx[tid];
    device const uchar* Q = queries + queryOffsets[readIdx];
    int qlen = int(queryLengths[readIdx]);
    int startPos = int(taskStartPos[tid]);
    int64_t minIntv = taskMinIntv[tid];
    short minSeedLen = minSeedLenBuf[0];
    int64_t sent = sentinel[0];

    uchar outIdx = 0;
    uint base_off = tid * MAX_OUT;

    // Check base at startPos
    int a = int(Q[startPos]);
    if (a >= 4) {
        outCount[tid] = 0;
        return;
    }

    // Initialize bidirectional SA interval for first character
    int64_t cur_k = bwtCounts[a];
    int64_t cur_l = bwtCounts[3 - a];
    int64_t cur_s = bwtCounts[a + 1] - bwtCounts[a];

    BidiInterval prevArr[MAX_PREV];
    int numPrev = 0;

    BidiInterval cur = {cur_k, cur_l, cur_s, startPos, startPos};

    // === Forward phase: extend rightward from startPos+1 ===
    for (int j = startPos + 1; j < qlen; j++) {
        int base = int(Q[j]);
        if (base >= 4) break;

        // Forward extension
        int64_t fk = cur.k, fl = cur.l, fs = cur.s;
        forward_ext(checkpoints, bwtCounts, sent, base, fk, fl, fs);

        // If interval size changed, record previous as candidate
        if (fs != cur.s && numPrev < MAX_PREV) {
            prevArr[numPrev] = cur;
            numPrev++;
        }

        // If new interval too small, stop
        if (fs < minIntv) {
            break;
        }

        cur = {fk, fl, fs, cur.m, j};
    }

    // Record final forward interval
    if (cur.s >= minIntv && numPrev < MAX_PREV) {
        prevArr[numPrev] = cur;
        numPrev++;
    }

    // Reverse prevArray (longest match first)
    for (int i = 0; i < numPrev / 2; i++) {
        BidiInterval tmp = prevArr[i];
        prevArr[i] = prevArr[numPrev - 1 - i];
        prevArr[numPrev - 1 - i] = tmp;
    }

    // === Backward phase: extend leftward from startPos-1 ===
    for (int bj = startPos - 1; bj >= 0; bj--) {
        int base = int(Q[bj]);
        if (base >= 4) break;

        int numCurr = 0;
        int64_t currS = -1;

        // First loop
        int p = 0;
        while (p < numPrev) {
            BidiInterval interval = prevArr[p];
            int64_t new_k, new_l, new_s;
            backward_ext(checkpoints, bwtCounts, sent,
                         interval.k, interval.l, interval.s,
                         base, new_k, new_l, new_s);

            // CONDITION 1: Can't extend further AND meets min length
            if (new_s < minIntv && (interval.n - interval.m + 1) >= int(minSeedLen)) {
                if (outIdx < MAX_OUT) {
                    outK[base_off + outIdx] = interval.k;
                    outL[base_off + outIdx] = interval.k + interval.s - 1;
                    outQBegin[base_off + outIdx] = short(interval.m);
                    outQEnd[base_off + outIdx] = short(interval.n + 1);
                    outIdx++;
                }
                break;
            }

            // CONDITION 2: Can extend AND hasn't seen this size yet
            if (new_s >= minIntv && new_s != currS) {
                currS = new_s;
                BidiInterval ni = {new_k, new_l, new_s, bj, interval.n};
                prevArr[numCurr] = ni;
                numCurr++;
                break;
            }
            p++;
        }

        // Second loop: continue from p+1
        p++;
        while (p < numPrev) {
            BidiInterval interval = prevArr[p];
            int64_t new_k, new_l, new_s;
            backward_ext(checkpoints, bwtCounts, sent,
                         interval.k, interval.l, interval.s,
                         base, new_k, new_l, new_s);

            if (new_s >= minIntv && new_s != currS) {
                currS = new_s;
                BidiInterval ni = {new_k, new_l, new_s, bj, interval.n};
                prevArr[numCurr] = ni;
                numCurr++;
            }
            p++;
        }

        numPrev = numCurr;
        if (numCurr == 0) break;
    }

    // Emit remaining intervals
    if (numPrev != 0) {
        BidiInterval interval = prevArr[0];
        if ((interval.n - interval.m + 1) >= int(minSeedLen) && outIdx < MAX_OUT) {
            outK[base_off + outIdx] = interval.k;
            outL[base_off + outIdx] = interval.k + interval.s - 1;
            outQBegin[base_off + outIdx] = short(interval.m);
            outQEnd[base_off + outIdx] = short(interval.n + 1);
            outIdx++;
        }
    }

    outCount[tid] = outIdx;
}
