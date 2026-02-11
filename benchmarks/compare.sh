#!/usr/bin/env bash
set -euo pipefail

# compare.sh — Field-by-field SAM comparison between bwa-mem2 and swiftbwa
# Usage: ./compare.sh [num_reads] [threads]
# Default: 10000 reads, 1 thread

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results"

NUM_READS="${1:-10000}"
THREADS="${2:-1}"
LABEL="${NUM_READS}_t${THREADS}"

# Tool paths
BWA_MEM2="/Users/jkane/Applications/bwa-mem2/bwa-mem2"
SWIFTBWA="$PROJECT_ROOT/.build/release/swiftbwa"
SAMTOOLS="/opt/homebrew/bin/samtools"

# Reference index
CHM13_INDEX="/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa.gz"

# Input FASTQs
R1="$DATA_DIR/hg002_${NUM_READS}_R1.fq"
R2="$DATA_DIR/hg002_${NUM_READS}_R2.fq"

# ─── Verify ──────────────────────────────────────────────────────────────────

echo "=== SAM Comparison: bwa-mem2 vs swiftbwa ==="
echo "Reads: $NUM_READS  Threads: $THREADS"
echo ""

for tool in "$BWA_MEM2" "$SWIFTBWA" "$SAMTOOLS"; do
    if [[ ! -x "$tool" ]]; then
        echo "ERROR: missing: $tool"
        exit 1
    fi
done
for f in "$R1" "$R2"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing FASTQ: $f (run setup.sh first)"
        exit 1
    fi
done

mkdir -p "$RESULTS_DIR"

# ─── Run aligners ───────────────────────────────────────────────────────────

TMP_BWA=$(mktemp /tmp/compare_bwa_XXXXXX.sam)
TMP_SWIFT=$(mktemp /tmp/compare_swift_XXXXXX.sam)
trap "rm -f $TMP_BWA $TMP_SWIFT" EXIT

echo "Running bwa-mem2..."
"$BWA_MEM2" mem -t "$THREADS" "$CHM13_INDEX" "$R1" "$R2" > "$TMP_BWA" 2>/dev/null

echo "Running swiftbwa..."
"$SWIFTBWA" mem -t "$THREADS" "$CHM13_INDEX" "$R1" "$R2" > "$TMP_SWIFT" 2>/dev/null

# ─── Extract primary alignments, sort by QNAME ──────────────────────────────

# Extract non-header, primary-only (FLAG & 0x900 == 0) lines, sort by QNAME
extract_primary() {
    local sam_file="$1"
    awk -F'\t' '!/^@/ && and($2, 0x900) == 0 { print }' "$sam_file" | sort -k1,1 -k2,2n
}

TMP_BWA_PRI=$(mktemp /tmp/compare_bwa_pri_XXXXXX.tsv)
TMP_SWIFT_PRI=$(mktemp /tmp/compare_swift_pri_XXXXXX.tsv)
trap "rm -f $TMP_BWA $TMP_SWIFT $TMP_BWA_PRI $TMP_SWIFT_PRI" EXIT

echo "Extracting primary alignments..."
extract_primary "$TMP_BWA" > "$TMP_BWA_PRI"
extract_primary "$TMP_SWIFT" > "$TMP_SWIFT_PRI"

BWA_PRI_COUNT=$(wc -l < "$TMP_BWA_PRI" | tr -d ' ')
SWIFT_PRI_COUNT=$(wc -l < "$TMP_SWIFT_PRI" | tr -d ' ')
echo "Primary records: bwa-mem2=$BWA_PRI_COUNT  swiftbwa=$SWIFT_PRI_COUNT"

# ─── Join and compare ────────────────────────────────────────────────────────

echo ""
echo "=== Comparison Results ==="
echo ""

TSV_OUT="$RESULTS_DIR/comparison_${LABEL}.tsv"

# The main comparison awk script:
# Reads both sorted files and joins on QNAME+read1/2 key, then compares fields.
awk -F'\t' '
BEGIN {
    OFS = "\t"
    # Counters
    bwa_total = 0; swift_total = 0
    bwa_mapped = 0; swift_mapped = 0
    both_mapped = 0
    bwa_only = 0; swift_only = 0
    pos_exact = 0; pos_5bp = 0; pos_diffchr = 0
    cigar_match = 0; cigar_total = 0
    mapq_exact = 0; mapq_total = 0; mapq_diff_sum = 0
    nm_match = 0; nm_total = 0
    flag_match = 0; flag_total = 0
    bwa_supp = 0; swift_supp = 0
    bwa_proper = 0; swift_proper = 0
}
# Read bwa-mem2 file (NR == FNR means first file)
NR == FNR {
    key = $1
    if (and($2, 0x40)) key = key "/1"
    else if (and($2, 0x80)) key = key "/2"
    bwa_flag[key] = $2
    bwa_rname[key] = $3
    bwa_pos[key] = $4
    bwa_mapq[key] = $5
    bwa_cigar[key] = $6
    bwa_tlen[key] = $9
    # Extract NM tag
    for (i = 12; i <= NF; i++) {
        if (substr($i, 1, 5) == "NM:i:") {
            bwa_nm[key] = substr($i, 6)
        }
    }
    bwa_total++
    if (and($2, 0x4) == 0) bwa_mapped++
    if (and($2, 0x2)) bwa_proper++
    next
}
# Read swiftbwa file (second file)
{
    key = $1
    if (and($2, 0x40)) key = key "/1"
    else if (and($2, 0x80)) key = key "/2"
    swift_total++
    if (and($2, 0x4) == 0) swift_mapped++
    if (and($2, 0x2)) swift_proper++

    # Extract NM
    swift_nm_val = ""
    for (i = 12; i <= NF; i++) {
        if (substr($i, 1, 5) == "NM:i:") {
            swift_nm_val = substr($i, 6)
        }
    }

    # Check if we have a bwa record for this key
    if (!(key in bwa_flag)) {
        swift_only_reads++
        next
    }

    bf = bwa_flag[key]
    sf = $2
    bwa_unmapped = and(bf, 0x4)
    swift_unmapped = and(sf, 0x4)

    if (bwa_unmapped && swift_unmapped) {
        # Both unmapped — nothing to compare
        next
    }
    if (!bwa_unmapped && swift_unmapped) {
        bwa_only++
        next
    }
    if (bwa_unmapped && !swift_unmapped) {
        swift_only++
        next
    }

    # Both mapped
    both_mapped++

    # Position concordance
    if (bwa_rname[key] == $3) {
        diff = bwa_pos[key] - $4
        if (diff < 0) diff = -diff
        if (diff == 0) pos_exact++
        else if (diff <= 5) pos_5bp++
    } else {
        pos_diffchr++
    }

    # CIGAR (only for same-position reads)
    if (bwa_rname[key] == $3 && bwa_pos[key] == $4) {
        cigar_total++
        if (bwa_cigar[key] == $6) cigar_match++
    }

    # MAPQ
    mapq_total++
    mdiff = bwa_mapq[key] - $5
    if (mdiff < 0) mdiff = -mdiff
    mapq_diff_sum += mdiff
    if (bwa_mapq[key] == $5) mapq_exact++

    # NM (only for same-position reads)
    if (bwa_rname[key] == $3 && bwa_pos[key] == $4) {
        if (bwa_nm[key] != "" && swift_nm_val != "") {
            nm_total++
            if (bwa_nm[key] == swift_nm_val) nm_match++
        }
    }

    # FLAG concordance (essential bits: 0x4|0x10|0x1|0x2|0x8|0x20|0x40|0x80)
    flag_total++
    mask = or(or(or(0x4, 0x10), or(0x1, 0x2)), or(or(0x8, 0x20), or(0x40, 0x80)))
    if (and(bf, mask) == and(sf, mask)) flag_match++
}
END {
    # Mapping rate
    bwa_map_pct = (bwa_total > 0) ? 100.0 * bwa_mapped / bwa_total : 0
    swift_map_pct = (swift_total > 0) ? 100.0 * swift_mapped / swift_total : 0
    printf "── Mapping Rate ──\n"
    printf "  bwa-mem2:  %d / %d (%.1f%%)\n", bwa_mapped, bwa_total, bwa_map_pct
    printf "  swiftbwa:  %d / %d (%.1f%%)\n", swift_mapped, swift_total, swift_map_pct
    printf "  Both mapped:   %d\n", both_mapped
    printf "  bwa-only:      %d\n", bwa_only
    printf "  swift-only:    %d\n", swift_only
    printf "\n"

    # Position concordance
    pos_compared = both_mapped
    pos_exact_pct = (pos_compared > 0) ? 100.0 * pos_exact / pos_compared : 0
    pos_5bp_pct = (pos_compared > 0) ? 100.0 * pos_5bp / pos_compared : 0
    pos_diff_pct = (pos_compared > 0) ? 100.0 * pos_diffchr / pos_compared : 0
    printf "── Position Concordance (both-mapped) ──\n"
    printf "  Exact match:    %d / %d (%.1f%%)\n", pos_exact, pos_compared, pos_exact_pct
    printf "  Within 5bp:     %d (%.1f%%)\n", pos_5bp, pos_5bp_pct
    printf "  Different chr:  %d (%.1f%%)\n", pos_diffchr, pos_diff_pct
    printf "\n"

    # CIGAR
    cigar_pct = (cigar_total > 0) ? 100.0 * cigar_match / cigar_total : 0
    printf "── CIGAR Agreement (same-position reads) ──\n"
    printf "  Exact match:  %d / %d (%.1f%%)\n", cigar_match, cigar_total, cigar_pct
    printf "\n"

    # MAPQ
    mapq_exact_pct = (mapq_total > 0) ? 100.0 * mapq_exact / mapq_total : 0
    mapq_mean_diff = (mapq_total > 0) ? mapq_diff_sum / mapq_total : 0
    printf "── MAPQ Agreement (both-mapped) ──\n"
    printf "  Exact match:    %d / %d (%.1f%%)\n", mapq_exact, mapq_total, mapq_exact_pct
    printf "  Mean abs diff:  %.2f\n", mapq_mean_diff
    printf "\n"

    # NM
    nm_pct = (nm_total > 0) ? 100.0 * nm_match / nm_total : 0
    printf "── NM Agreement (same-position reads) ──\n"
    printf "  Exact match:  %d / %d (%.1f%%)\n", nm_match, nm_total, nm_pct
    printf "\n"

    # FLAG
    flag_pct = (flag_total > 0) ? 100.0 * flag_match / flag_total : 0
    printf "── FLAG Concordance (essential bits) ──\n"
    printf "  Match:  %d / %d (%.1f%%)\n", flag_match, flag_total, flag_pct
    printf "\n"

    # Proper pair
    bwa_pp_pct = (bwa_total > 0) ? 100.0 * bwa_proper / bwa_total : 0
    swift_pp_pct = (swift_total > 0) ? 100.0 * swift_proper / swift_total : 0
    printf "── Proper Pair Rate ──\n"
    printf "  bwa-mem2:  %.1f%%\n", bwa_pp_pct
    printf "  swiftbwa:  %.1f%%\n", swift_pp_pct
    printf "\n"
}
' "$TMP_BWA_PRI" "$TMP_SWIFT_PRI" | tee "$TSV_OUT"

# Supplementary counts from full SAM
BWA_SUPP=$(awk -F'\t' '!/^@/ && and($2, 0x800) { c++ } END { print c+0 }' "$TMP_BWA")
SWIFT_SUPP=$(awk -F'\t' '!/^@/ && and($2, 0x800) { c++ } END { print c+0 }' "$TMP_SWIFT")

echo "── Supplementary Alignments ──"
echo "  bwa-mem2:  $BWA_SUPP"
echo "  swiftbwa:  $SWIFT_SUPP"
echo ""

echo "Results saved to: $TSV_OUT"
