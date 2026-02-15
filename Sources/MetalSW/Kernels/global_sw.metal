#include <metal_stdlib>
using namespace metal;

/// Banded global (Needleman-Wunsch) alignment + traceback kernel.
/// One GPU thread per alignment task. Produces packed CIGAR.
///
/// Matches GlobalAligner.align() + traceback() from GlobalAligner.swift.
/// Direction bit encoding (same as Swift):
///   bit 0: E extend (0=opened from H, 1=extended)
///   bit 1: F extend (0=opened from H, 1=extended)
///   bit 2: H came from E
///   bit 3: H came from F

constant int MAX_CIGAR_OPS = 512;

kernel void global_sw(
    device const uchar*  queries       [[buffer(0)]],   // contiguous packed queries
    device const uint*   queryOffsets  [[buffer(1)]],   // byte offset per task
    device const ushort* queryLengths  [[buffer(2)]],   // length per task
    device const uchar*  targets       [[buffer(3)]],   // contiguous packed targets
    device const uint*   targetOffsets [[buffer(4)]],   // byte offset per task
    device const ushort* targetLengths [[buffer(5)]],   // length per task
    device const char*   scoringMatrix [[buffer(6)]],   // 5x5 flat Int8
    device const short*  params        [[buffer(7)]],   // [oIns,eIns,oDel,eDel]
    device const short*  trueScores    [[buffer(8)]],   // expected score per task
    device const short*  initialW      [[buffer(9)]],   // initial bandwidth per task
    device       int*    cigarOut      [[buffer(10)]],  // packed CIGAR per task
    device       int*    cigarLengths  [[buffer(11)]],  // # CIGAR ops per task
    device       int*    scores        [[buffer(12)]],  // output score per task
    device       uchar*  workspace     [[buffer(13)]],  // per-task workspace
    device const uint*   wsOffsets     [[buffer(14)]],  // workspace byte offset per task
    uint tid [[thread_position_in_grid]]
)
{
    const int qlen = int(queryLengths[tid]);
    const int tlen = int(targetLengths[tid]);

    // Output pointers
    device int* mycigar = cigarOut + tid * MAX_CIGAR_OPS;

    if (qlen == 0 || tlen == 0) {
        cigarLengths[tid] = 0;
        scores[tid] = 0;
        return;
    }

    // Scoring parameters
    const int oIns = int(params[0]);   // gap open insertion
    const int eIns = int(params[1]);   // gap extend insertion
    const int oDel = int(params[2]);   // gap open deletion
    const int eDel = int(params[3]);   // gap extend deletion
    const int oeIns = oIns + eIns;
    const int oeDel = oDel + eDel;
    const int negInf = -1073741824;    // INT_MIN/2 equivalent, avoids overflow

    device const uchar* query  = queries + queryOffsets[tid];
    device const uchar* target = targets + targetOffsets[tid];

    const int trueScore = int(trueScores[tid]);
    int w = int(initialW[tid]);

    // Ensure bandwidth covers length difference
    int lenDiff = qlen - tlen;
    if (lenDiff < 0) lenDiff = -lenDiff;
    if (w < lenDiff) w = lenDiff;

    // Workspace layout:
    //   qp:  5 * qlen * sizeof(int)     = 20 * qlen bytes
    //   H:   (qlen+1) * sizeof(int)     = 4 * (qlen+1) bytes
    //   E:   (qlen+1) * sizeof(int)     = 4 * (qlen+1) bytes
    //   z:   tlen * nCol * sizeof(uchar) = tlen * nCol bytes (variable with w)
    // We allocate max z for worst-case w (after retries w can grow up to 8*initialW).
    device uchar* ws = workspace + wsOffsets[tid];
    device int* qp = (device int*)ws;
    device int* H  = qp + 5 * qlen;
    device int* E  = H + (qlen + 1);
    device uchar* zBase = (device uchar*)(E + (qlen + 1));

    // Build query profile: qp[k * qlen + j] = scoringMatrix[k * 5 + query[j]]
    for (int k = 0; k < 5; k++) {
        for (int j = 0; j < qlen; j++) {
            qp[k * qlen + j] = int(scoringMatrix[k * 5 + int(query[j])]);
        }
    }

    // Retry loop: double bandwidth if score doesn't match trueScore
    int finalScore = negInf;
    int retries = 0;
    const int maxRetries = 3;

    while (true) {
        int nCol = qlen < (2 * w + 1) ? qlen : (2 * w + 1);
        device uchar* z = zBase;

        // Initialize DP arrays
        for (int j = 0; j <= qlen; j++) {
            H[j] = negInf;
            E[j] = negInf;
        }
        H[0] = 0;
        int initEnd = w < qlen ? w : qlen;
        for (int j = 1; j <= initEnd; j++) {
            H[j] = -(oeIns + eIns * (j - 1));
        }

        // Initialize z
        for (int idx = 0; idx < tlen * nCol; idx++) {
            z[idx] = 0;
        }

        // Main DP loop — banded NW matching GlobalAligner.swift lines 89-141
        for (int i = 0; i < tlen; i++) {
            int ii = int(i);
            int beg = ii - w;
            if (beg < 0) beg = 0;
            int end = ii + w + 1;
            if (end > qlen) end = qlen;

            int targetBase = min(int(target[i]), 4);
            int qpRow = targetBase * qlen;

            int f = negInf;
            int hDiag = H[beg];

            // Left boundary
            if (beg == 0) {
                H[beg] = -(oeDel + eDel * ii);
            } else {
                H[beg] = negInf;
            }

            device uchar* zi = z + (i * nCol - beg);

            for (int j = beg; j < end; j++) {
                int j1 = j + 1;

                // Diagonal: match/mismatch
                int hCur = hDiag + qp[qpRow + j];

                // Save H(i-1,j) for next iteration's diagonal
                hDiag = H[j1];

                // E: deletion (gap in query, vertical move)
                int eOpen = hDiag - oeDel;
                int eExtend = E[j1] - eDel;
                uchar eFlag = (eOpen > eExtend) ? 0 : 0x01;
                int eVal = (eOpen > eExtend) ? eOpen : eExtend;
                E[j1] = eVal;

                // F: insertion (gap in target, horizontal move)
                int fOpen = H[j] - oeIns;
                int fExtend = f - eIns;
                uchar fFlag = (fOpen > fExtend) ? 0 : 0x02;
                int fVal = (fOpen > fExtend) ? fOpen : fExtend;
                f = fVal;

                // H = max(hCur, E, F)
                uchar fromE = (eVal > hCur) ? 0x04 : 0;
                hCur = (eVal > hCur) ? eVal : hCur;
                uchar fromF = (fVal > hCur) ? 0x08 : 0;
                hCur = (fVal > hCur) ? fVal : hCur;

                H[j1] = hCur;
                zi[j] = eFlag | fFlag | fromE | fromF;
            }
        }

        finalScore = H[qlen];

        if (finalScore >= trueScore || retries >= maxRetries) {
            // Traceback
            // Build CIGAR in reverse into a local stack, then reverse into output
            int tmpCigar[MAX_CIGAR_OPS];  // thread-local stack
            int tmpCount = 0;

            int nCol_tb = qlen < (2 * w + 1) ? qlen : (2 * w + 1);
            device uchar* z_tb = zBase;

            int ci = tlen - 1;  // target position
            int cj = qlen - 1;  // query position
            int state = 0;      // 0=M, 1=E(deletion), 2=F(insertion)

            while ((ci >= 0 || cj >= 0) && tmpCount < MAX_CIGAR_OPS) {
                if (ci < 0) {
                    // Only query bases remain -> insertions
                    // Pack: op=1 (insertion), length=1
                    if (tmpCount > 0 && (tmpCount > 0) && ((tmpCigar[tmpCount - 1] & 0xF) == 1)) {
                        tmpCigar[tmpCount - 1] += (1 << 4);
                    } else {
                        tmpCigar[tmpCount++] = (1 << 4) | 1;
                    }
                    cj--;
                    continue;
                }
                if (cj < 0) {
                    // Only target bases remain -> deletions
                    if (tmpCount > 0 && ((tmpCigar[tmpCount - 1] & 0xF) == 2)) {
                        tmpCigar[tmpCount - 1] += (1 << 4);
                    } else {
                        tmpCigar[tmpCount++] = (1 << 4) | 2;
                    }
                    ci--;
                    continue;
                }

                int beg_tb = ci - w;
                if (beg_tb < 0) beg_tb = 0;
                int col = cj - beg_tb;

                if (col < 0 || col >= nCol_tb) {
                    break;
                }

                uchar d = z_tb[ci * nCol_tb + col];

                if (state == 0) {
                    // M state
                    int hDir = (int(d) >> 2) & 0x03;
                    if (hDir == 0) {
                        // Diagonal: emit M
                        if (tmpCount > 0 && ((tmpCigar[tmpCount - 1] & 0xF) == 0)) {
                            tmpCigar[tmpCount - 1] += (1 << 4);
                        } else {
                            tmpCigar[tmpCount++] = (1 << 4) | 0;
                        }
                        ci--;
                        cj--;
                    } else if (hDir == 1) {
                        state = 1;  // switch to E (deletion)
                    } else {
                        state = 2;  // switch to F (insertion)
                    }
                } else if (state == 1) {
                    // E state: deletion (gap in query, consume target)
                    if (tmpCount > 0 && ((tmpCigar[tmpCount - 1] & 0xF) == 2)) {
                        tmpCigar[tmpCount - 1] += (1 << 4);
                    } else {
                        tmpCigar[tmpCount++] = (1 << 4) | 2;
                    }
                    if ((d & 0x01) == 0) {
                        state = 0;  // opened from H
                    }
                    ci--;
                } else {
                    // F state: insertion (gap in target, consume query)
                    if (tmpCount > 0 && ((tmpCigar[tmpCount - 1] & 0xF) == 1)) {
                        tmpCigar[tmpCount - 1] += (1 << 4);
                    } else {
                        tmpCigar[tmpCount++] = (1 << 4) | 1;
                    }
                    if ((d & 0x02) == 0) {
                        state = 0;  // opened from H
                    }
                    cj--;
                }
            }

            // Reverse CIGAR into output
            int outCount = tmpCount;
            for (int r = 0; r < outCount; r++) {
                mycigar[r] = tmpCigar[outCount - 1 - r];
            }
            cigarLengths[tid] = outCount;
            scores[tid] = finalScore;
            return;
        }

        // Retry with doubled bandwidth
        w *= 2;
        retries++;
    }
}
