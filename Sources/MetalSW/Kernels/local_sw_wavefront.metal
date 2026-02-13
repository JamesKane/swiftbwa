#include <metal_stdlib>
using namespace metal;

/// Local Smith-Waterman using anti-diagonal wavefront parallelism.
/// One SIMD group (32 threads) cooperates on each alignment task.
/// Anti-diagonals are processed sequentially; within each anti-diagonal,
/// 32 query positions are computed in parallel using simd_shuffle.
/// For queries > 32bp, tiles of 32 positions are processed sequentially
/// with inter-tile boundary values passed via simd_broadcast.
///
/// Results layout per task (3 ints):
///   [0] score  (-1 if no alignment found)
///   [1] queryEnd
///   [2] targetEnd

#define SIMD_W 32
#define MAX_TILES 5   // supports up to 160bp queries

kernel void local_sw_wavefront(
    device const uchar*  queries       [[buffer(0)]],
    device const uint*   queryOffsets  [[buffer(1)]],
    device const ushort* queryLengths  [[buffer(2)]],
    device const uchar*  targets       [[buffer(3)]],
    device const uint*   targetOffsets [[buffer(4)]],
    device const ushort* targetLengths [[buffer(5)]],
    device const char*   scoringMatrix [[buffer(6)]],  // 5x5 signed
    device const short*  params        [[buffer(7)]],  // [gapOE, gapE]
    device       int*    results       [[buffer(8)]],
    device const uint&   numTasks      [[buffer(9)]],
    uint lane    [[thread_index_in_simdgroup]],
    uint sg_idx  [[simdgroup_index_in_threadgroup]],
    uint tg_idx  [[threadgroup_position_in_grid]],
    uint sg_per_tg [[simdgroups_per_threadgroup]]
) {
    uint task = tg_idx * sg_per_tg + sg_idx;
    if (task >= numTasks) return;

    uint qlen = queryLengths[task];
    uint tlen = targetLengths[task];
    device const uchar* Q = queries + queryOffsets[task];
    device const uchar* T = targets + targetOffsets[task];

    if (qlen == 0 || tlen == 0) {
        if (lane == 0) {
            results[task * 3]     = -1;
            results[task * 3 + 1] = 0;
            results[task * 3 + 2] = 0;
        }
        return;
    }

    int gapOE = params[0];  // gap open + extend
    int gapE  = params[1];  // gap extend only

    uint n_tiles = (qlen + SIMD_W - 1) / SIMD_W;
    if (n_tiles > MAX_TILES) n_tiles = MAX_TILES;

    // Per-tile DP state in registers
    // hp[t]  = H value for my query position in tile t, from anti-diag d-1
    // hp2[t] = H value from anti-diag d-2
    // ep[t]  = E value from anti-diag d-1
    // fp[t]  = F value from anti-diag d-1
    int hp[MAX_TILES]  = {0, 0, 0, 0, 0};
    int hp2[MAX_TILES] = {0, 0, 0, 0, 0};
    int ep[MAX_TILES]  = {0, 0, 0, 0, 0};
    int fp[MAX_TILES]  = {0, 0, 0, 0, 0};

    int best = 0, best_qi = -1, best_tj = -1;

    uint n_adiag = qlen + tlen - 1;

    for (uint d = 0; d < n_adiag; d++) {
        // Save thread-31 boundary values BEFORE updating any tile on this anti-diagonal.
        // Tile t+1's thread 0 needs tile t's thread 31 values from BEFORE this anti-diagonal.
        int bh2[MAX_TILES], bh1[MAX_TILES], bf[MAX_TILES];
        for (uint t = 0; t < n_tiles; t++) {
            bh2[t] = simd_broadcast(hp2[t], 31);
            bh1[t] = simd_broadcast(hp[t], 31);
            bf[t]  = simd_broadcast(fp[t], 31);
        }

        for (uint t = 0; t < n_tiles; t++) {
            uint qi = t * SIMD_W + lane;
            int  tj = int(d) - int(qi);

            // Load this tile's state from previous anti-diagonal
            int h_old = hp2[t];
            int h_cur = hp[t];
            int e_cur = ep[t];
            int f_cur = fp[t];

            // Shuffle: get values from thread lane-1 (query position qi-1 within this tile)
            // H diagonal: H[qi-1] from 2 anti-diags ago (for H[i-1][j-1])
            int h_diag = simd_shuffle_up(h_old, 1);
            // H upper: H[qi-1] from 1 anti-diag ago (for F[i-1][j])
            int h_up   = simd_shuffle_up(h_cur, 1);
            // F upper: F[qi-1] from 1 anti-diag ago (for F[i-1][j])
            int f_up   = simd_shuffle_up(f_cur, 1);

            // Thread 0 gets boundary from previous tile (or zero for tile 0)
            if (lane == 0) {
                if (t > 0) {
                    h_diag = bh2[t - 1];
                    h_up   = bh1[t - 1];
                    f_up   = bf[t - 1];
                } else {
                    h_diag = 0;
                    h_up   = 0;
                    f_up   = 0;
                }
            }

            int h_new, e_new, f_new;

            if (qi < qlen && tj >= 0 && uint(tj) < tlen) {
                int score = int(scoringMatrix[Q[qi] * 5 + T[tj]]);

                // E: horizontal gap extension (gap in query, advancing target)
                // E[i][j] = max(H[i][j-1] - gapOE, E[i][j-1] - gapE)
                // On anti-diag d, "j-1" is from anti-diag d-1, same thread
                e_new = max(h_cur - gapOE, e_cur - gapE);
                if (e_new < 0) e_new = 0;

                // F: vertical gap extension (gap in target, advancing query)
                // F[i][j] = max(H[i-1][j] - gapOE, F[i-1][j] - gapE)
                // On anti-diag d, "i-1" is from anti-diag d-1, thread lane-1
                f_new = max(h_up - gapOE, f_up - gapE);
                if (f_new < 0) f_new = 0;

                // H: best of diagonal + match, E, F, or 0 (local alignment)
                h_new = max(max(h_diag + score, e_new), max(f_new, 0));

                if (h_new > best) {
                    best = h_new;
                    best_qi = int(qi);
                    best_tj = tj;
                }
            } else {
                h_new = 0;
                e_new = 0;
                f_new = 0;
            }

            // Store updated state for next anti-diagonal
            hp2[t] = h_cur;    // d-1 becomes d-2
            hp[t]  = h_new;    // new H becomes d-1
            ep[t]  = e_new;
            fp[t]  = f_new;
        }
    }

    // Reduce max score across SIMD group
    int group_best = simd_max(best);

    // Find lowest lane holding the max (deterministic tiebreak)
    uint winner_rank = (best == group_best && best > 0) ? lane : SIMD_W;
    uint winner = simd_min(winner_rank);

    int final_qi = simd_broadcast(best_qi, ushort(winner));
    int final_tj = simd_broadcast(best_tj, ushort(winner));

    if (lane == 0) {
        device int* out = results + task * 3;
        if (group_best <= 0) {
            out[0] = -1;
            out[1] = 0;
            out[2] = 0;
        } else {
            out[0] = group_best;
            out[1] = final_qi;
            out[2] = final_tj;
        }
    }
}
