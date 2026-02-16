#include <metal_stdlib>
using namespace metal;

/// 16-bit banded Smith-Waterman kernel — scalar per-element DP.
/// One GPU thread per alignment task (inter-sequence parallelism).
/// Uses signed 16-bit (no bias needed). Handles overflow cases from 8-bit kernel.
///
/// Results layout per task (6 ints):
///   [0] score
///   [1] queryEnd
///   [2] targetEnd
///   [3] globalTargetEnd
///   [4] globalScore
///   [5] maxOff
kernel void banded_sw16(
    device const uchar*   queries        [[buffer(0)]],
    device const uint*    queryOffsets    [[buffer(1)]],
    device const ushort*  queryLengths   [[buffer(2)]],
    device const uchar*   targets        [[buffer(3)]],
    device const uint*    targetOffsets   [[buffer(4)]],
    device const ushort*  targetLengths  [[buffer(5)]],
    device const char*    scoringMatrix  [[buffer(6)]],  // 5x5 flat
    device const short*   params         [[buffer(7)]],  // [oIns,eIns,oDel,eDel,zDrop,_,w]
    device const short*   h0Values       [[buffer(8)]],
    device       int*     results        [[buffer(9)]],
    device       uchar*   workspace      [[buffer(10)]],
    device const uint*    wsOffsets      [[buffer(11)]],
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
    const short w     = params[6];
    const short h0    = h0Values[tid];

    device const uchar* query  = queries + queryOffsets[tid];
    device const uchar* target = targets + targetOffsets[tid];

    // Workspace: H[qlen] + E[qlen] + profile[5*qlen], each element is short (2 bytes)
    device short* H = (device short*)(workspace + wsOffsets[tid]);
    device short* E = H + qlen;
    device short* profile = E + qlen;

    // Build query profile (signed, no bias)
    for (int k = 0; k < 5; k++) {
        for (int j = 0; j < qlen; j++) {
            profile[k * qlen + j] = short(scoringMatrix[k * 5 + int(query[j])]);
        }
    }

    // Initialize H
    for (int j = 0; j < qlen; j++) {
        if (j + 1 <= int(w)) {
            int val = int(h0) - int(oIns) - int(eIns) * (j + 1);
            H[j] = val > 0 ? short(min(val, 32767)) : 0;
        } else {
            H[j] = 0;
        }
        E[j] = 0;
    }

    int maxScore = int(h0);
    int maxI = -1;
    int maxJ = -1;
    int gScore = -1;
    int gTle = -1;
    int maxOff = 0;

    for (int i = 0; i < tlen; i++) {
        int targetBase = min(int(target[i]), 4);
        device short* prof = profile + targetBase * qlen;

        short hDiag = (i == 0) ? short(clamp(max(0, int(h0)), 0, 32767)) : 0;
        short f = 0;
        int rowMax = 0;
        int rowMaxJ = -1;

        for (int j = 0; j < qlen; j++) {
            // Diagonal + profile
            // Zero out when hDiag == 0 to prevent alignment restart (bwa-mem2: M = M? M + q[j] : 0)
            int hNew = (hDiag != 0) ? int(hDiag) + int(prof[j]) : 0;
            if (hNew < 0) hNew = 0;

            short hPrev = H[j];
            hDiag = hPrev;

            // E: insertion
            int eOpen = int(hPrev) - int(oIns);
            int eVal = max(eOpen, int(E[j])) - int(eIns);
            if (eVal < 0) eVal = 0;
            E[j] = short(clamp(eVal, 0, 32767));

            hNew = max(hNew, int(E[j]));

            // F: deletion
            hNew = max(hNew, int(f));

            H[j] = short(clamp(hNew, 0, 32767));

            // Update F for next j
            int fOpen = hNew - int(oDel);
            int fVal = max(fOpen, int(f)) - int(eDel);
            if (fVal < 0) fVal = 0;
            f = short(clamp(fVal, 0, 32767));

            if (hNew > rowMax) { rowMax = hNew; rowMaxJ = j; }
        }

        // Early termination when row max is 0 (bwa-mem2 ksw.cpp:511)
        if (rowMax == 0) break;

        // Global score
        int hAtEnd = int(H[qlen - 1]);
        if (hAtEnd > gScore) {
            gScore = hAtEnd;
            gTle = i;
        }

        // Track overall maximum and z-drop (bwa-mem2 ksw.cpp:512-520)
        if (rowMax > maxScore) {
            maxScore = rowMax;
            maxI = i;
            maxJ = rowMaxJ;
            int off = abs(rowMaxJ - i);
            if (off > maxOff) maxOff = off;
        } else if (int(zDrop) > 0) {
            int deltaI = i - maxI;
            int deltaJ = rowMaxJ - maxJ;
            if (deltaI > deltaJ) {
                if (maxScore - rowMax - (deltaI - deltaJ) * int(eDel) > int(zDrop)) break;
            } else {
                if (maxScore - rowMax - (deltaJ - deltaI) * int(eIns) > int(zDrop)) break;
            }
        }
    }

    out[0] = maxScore;
    out[1] = maxJ + 1;
    out[2] = maxI + 1;
    out[3] = gTle + 1;
    out[4] = gScore;
    out[5] = maxOff;
}
