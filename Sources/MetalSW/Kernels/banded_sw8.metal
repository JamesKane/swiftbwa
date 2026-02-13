#include <metal_stdlib>
using namespace metal;

/// 8-bit banded Smith-Waterman kernel — scalar per-element DP.
/// One GPU thread per alignment task (inter-sequence parallelism).
/// Uses unsigned 8-bit with bias to keep values positive.
///
/// Results layout per task (6 ints):
///   [0] score  (INT_MIN = overflow, caller falls back to 16-bit)
///   [1] queryEnd
///   [2] targetEnd
///   [3] globalTargetEnd
///   [4] globalScore
///   [5] maxOff
kernel void banded_sw8(
    device const uchar*   queries        [[buffer(0)]],
    device const uint*    queryOffsets    [[buffer(1)]],
    device const ushort*  queryLengths   [[buffer(2)]],
    device const uchar*   targets        [[buffer(3)]],
    device const uint*    targetOffsets   [[buffer(4)]],
    device const ushort*  targetLengths  [[buffer(5)]],
    device const char*    scoringMatrix  [[buffer(6)]],  // 5x5 flat
    device const short*   params         [[buffer(7)]],  // [oIns,eIns,oDel,eDel,zDrop,bias,w]
    device const short*   h0Values       [[buffer(8)]],
    device       int*     results        [[buffer(9)]],
    device       uchar*   workspace      [[buffer(10)]],
    device const uint*    wsOffsets      [[buffer(11)]],  // workspace byte offset per task
    uint tid [[thread_position_in_grid]]
)
{
    const int qlen = int(queryLengths[tid]);
    const int tlen = int(targetLengths[tid]);
    device int* out = results + tid * 6;

    if (qlen == 0 || tlen == 0) {
        out[0] = 0; out[1] = 0; out[2] = 0;
        out[3] = 0; out[4] = 0; out[5] = 0;
        return;
    }

    const short oIns  = params[0];
    const short eIns  = params[1];
    const short oDel  = params[2];
    const short eDel  = params[3];
    const short zDrop = params[4];
    const short bias  = params[5];
    const short w     = params[6];
    const short h0    = h0Values[tid];

    device const uchar* query  = queries + queryOffsets[tid];
    device const uchar* target = targets + targetOffsets[tid];

    // Workspace: H[qlen] + E[qlen] as uint8 arrays
    device uchar* H = workspace + wsOffsets[tid];
    device uchar* E = H + qlen;

    // Build query profile with bias (5 * qlen bytes, after H and E)
    device uchar* profile = E + qlen;

    for (int k = 0; k < 5; k++) {
        for (int j = 0; j < qlen; j++) {
            char sc = scoringMatrix[k * 5 + int(query[j])];
            int val = int(sc) + int(bias);
            profile[k * qlen + j] = uchar(clamp(val, 0, 255));
        }
    }

    // Initialize H
    for (int j = 0; j < qlen; j++) {
        if (j + 1 <= int(w)) {
            int val = int(h0) - int(oIns) - int(eIns) * (j + 1);
            H[j] = val > 0 ? uchar(min(val, 255)) : 0;
        } else {
            H[j] = 0;
        }
        E[j] = 0;
    }

    uchar maxScore = uchar(clamp(int(h0), 0, 255));
    int maxI = -1;
    int maxJ = -1;
    int maxIE = -1;
    int gScore = -1;
    int gTle = -1;
    int maxOff = 0;

    for (int i = 0; i < tlen; i++) {
        int targetBase = min(int(target[i]), 4);
        device uchar* prof = profile + targetBase * qlen;

        // Process each query position
        uchar hDiag = (i == 0) ? uchar(clamp(max(0, int(h0)), 0, 255)) : 0;
        uchar f = 0;
        uchar rowMax = 0;

        for (int j = 0; j < qlen; j++) {
            // Diagonal + profile - bias
            int hNew = int(hDiag) + int(prof[j]) - int(bias);
            if (hNew < 0) hNew = 0;

            uchar hPrev = H[j];
            hDiag = hPrev;  // save for next j's diagonal

            // E: insertion (gap in target)
            int eOpen = max(int(hPrev), int(oIns)) - int(oIns);
            int eVal = max(eOpen, int(E[j]));
            eVal = max(eVal, int(eIns)) - int(eIns);
            E[j] = uchar(clamp(eVal, 0, 255));

            hNew = max(hNew, int(E[j]));

            // F: deletion (gap in query), carried from previous j
            hNew = max(hNew, int(f));

            if (hNew > 250) {
                // Overflow — signal caller to use 16-bit
                out[0] = -2147483647;
                out[1] = 0; out[2] = 0; out[3] = 0; out[4] = 0; out[5] = 0;
                return;
            }

            H[j] = uchar(hNew);

            // Update F for next j
            int fOpen = max(hNew, int(oDel)) - int(oDel);
            f = uchar(max(max(fOpen, int(f)), int(eDel)) - int(eDel));

            if (uchar(hNew) > rowMax) rowMax = uchar(hNew);

            if (uchar(hNew) > maxScore) {
                maxScore = uchar(hNew);
                maxI = i;
                maxJ = j;
                maxIE = i;
                int off = abs(maxI - maxJ);
                if (off > maxOff) maxOff = off;
            }
        }

        // Global score (alignment consuming entire query)
        int hAtEnd = int(H[qlen - 1]);
        if (hAtEnd > gScore) {
            gScore = hAtEnd;
            gTle = i;
        }

        // Z-dropoff
        if (int(maxScore) - int(rowMax) > int(zDrop) && i - maxIE > int(w)) {
            break;
        }
    }

    out[0] = int(maxScore);
    out[1] = maxJ + 1;
    out[2] = maxI + 1;
    out[3] = gTle + 1;
    out[4] = gScore;
    out[5] = maxOff;
}
