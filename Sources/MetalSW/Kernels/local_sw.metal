#include <metal_stdlib>
using namespace metal;

/// Local Smith-Waterman kernel — single-pass (forward or reverse).
/// One GPU thread per alignment task (inter-sequence parallelism).
/// Uses unsigned 8-bit with bias. Returns (score, qEnd, tEnd).
/// The two-pass start position recovery is handled by the Swift dispatcher
/// (forward batch → reverse batch).
///
/// Results layout per task (3 ints):
///   [0] score  (-1 if overflow)
///   [1] queryEnd
///   [2] targetEnd
kernel void local_sw(
    device const uchar*   queries        [[buffer(0)]],
    device const uint*    queryOffsets    [[buffer(1)]],
    device const ushort*  queryLengths   [[buffer(2)]],
    device const uchar*   targets        [[buffer(3)]],
    device const uint*    targetOffsets   [[buffer(4)]],
    device const ushort*  targetLengths  [[buffer(5)]],
    device const char*    scoringMatrix  [[buffer(6)]],  // 5x5 flat
    device const short*   params         [[buffer(7)]],  // [gapOE, gapE, bias]
    device       int*     results        [[buffer(8)]],
    device       uchar*   workspace      [[buffer(9)]],
    device const uint*    wsOffsets      [[buffer(10)]],
    uint tid [[thread_position_in_grid]]
)
{
    const int qlen = int(queryLengths[tid]);
    const int tlen = int(targetLengths[tid]);
    device int* out = results + tid * 3;

    if (qlen == 0 || tlen == 0) {
        out[0] = 0; out[1] = -1; out[2] = -1;
        return;
    }

    const short gapOE = params[0];
    const short gapE  = params[1];
    const uchar bias  = uchar(params[2]);

    device const uchar* query  = queries + queryOffsets[tid];
    device const uchar* target = targets + targetOffsets[tid];

    // Workspace: H[qlen] + E[qlen] + profile[5*qlen]
    device uchar* H = workspace + wsOffsets[tid];
    device uchar* E = H + qlen;
    device uchar* profile = E + qlen;

    // Build query profile with bias
    for (int k = 0; k < 5; k++) {
        for (int j = 0; j < qlen; j++) {
            char sc = scoringMatrix[k * 5 + int(query[j])];
            int val = int(sc) + int(bias);
            profile[k * qlen + j] = uchar(clamp(val, 0, 255));
        }
    }

    // Initialize H, E to zero
    for (int j = 0; j < qlen; j++) {
        H[j] = 0;
        E[j] = 0;
    }

    uchar maxScore = 0;
    int maxI = -1;
    int maxJ = -1;

    for (int i = 0; i < tlen; i++) {
        int targetBase = min(int(target[i]), 4);
        device uchar* prof = profile + targetBase * qlen;

        uchar hDiag = 0;
        uchar f = 0;

        for (int j = 0; j < qlen; j++) {
            // Diagonal + profile - bias
            int hNew = int(hDiag) + int(prof[j]) - int(bias);
            if (hNew < 0) hNew = 0;

            uchar hPrev = H[j];
            hDiag = hPrev;

            // E: gap in query (vertical)
            int eOpen = max(int(hPrev), int(gapOE)) - int(gapOE);
            int eVal = max(eOpen, int(E[j]));
            eVal = max(eVal, int(gapE)) - int(gapE);
            E[j] = uchar(clamp(eVal, 0, 255));

            hNew = max(hNew, int(E[j]));
            hNew = max(hNew, int(f));

            if (hNew > 250) {
                // Overflow
                out[0] = -1;
                out[1] = -1;
                out[2] = -1;
                return;
            }

            H[j] = uchar(hNew);

            // F: gap in target (horizontal)
            int fOpen = max(hNew, int(gapOE)) - int(gapOE);
            f = uchar(max(max(fOpen, int(f)), int(gapE)) - int(gapE));

            if (uchar(hNew) > maxScore) {
                maxScore = uchar(hNew);
                maxI = i;
                maxJ = j;
            }
        }
    }

    out[0] = int(maxScore);
    out[1] = maxJ;
    out[2] = maxI;
}
