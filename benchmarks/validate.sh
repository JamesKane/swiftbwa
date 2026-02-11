#!/usr/bin/env bash
set -euo pipefail

# validate.sh — GATK-based validation of SAM output correctness
# Usage: ./validate.sh [num_reads]
# Default: 10000 reads

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results"

NUM_READS="${1:-10000}"

# Tool paths
BWA_MEM2="/Users/jkane/Applications/bwa-mem2/bwa-mem2"
SWIFTBWA="$PROJECT_ROOT/.build/release/swiftbwa"
SAMTOOLS="/opt/homebrew/bin/samtools"
GATK="/Users/jkane/Applications/gatk-4.6.2.0/gatk"

# Reference (GATK needs uncompressed + .fai + .dict)
CHM13_REF="/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa"
CHM13_INDEX="/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa.gz"

# Input FASTQs
R1="$DATA_DIR/hg002_${NUM_READS}_R1.fq"
R2="$DATA_DIR/hg002_${NUM_READS}_R2.fq"

# ─── Verify ──────────────────────────────────────────────────────────────────

echo "=== GATK Validation: bwa-mem2 vs swiftbwa ==="
echo "Reads: $NUM_READS"
echo ""

fail=0
for tool in "$BWA_MEM2" "$SWIFTBWA" "$SAMTOOLS" "$GATK"; do
    if [[ ! -x "$tool" && ! -f "$tool" ]]; then
        echo "ERROR: missing: $tool"
        fail=1
    fi
done

# GATK needs uncompressed reference with .fai and .dict
for f in "$CHM13_REF" "$CHM13_REF.fai"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing reference file: $f"
        fail=1
    fi
done

# Check for .dict (could be chm13v2.0.dict or chm13v2.0.fa.dict)
DICT_FILE="/Users/jkane/Genomics/chm13v2.0/chm13v2.0.dict"
if [[ ! -f "$DICT_FILE" ]]; then
    DICT_FILE="$CHM13_REF.dict"
    if [[ ! -f "$DICT_FILE" ]]; then
        echo "ERROR: missing sequence dictionary (.dict)"
        fail=1
    fi
fi

for f in "$R1" "$R2"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing FASTQ: $f (run setup.sh first)"
        fail=1
    fi
done

if [[ $fail -ne 0 ]]; then
    echo "Fix errors above before running validation."
    exit 1
fi

echo "All tools and data verified."
mkdir -p "$RESULTS_DIR"

# ─── Run aligners ───────────────────────────────────────────────────────────

TMP_DIR=$(mktemp -d /tmp/validate_XXXXXX)
trap "rm -rf $TMP_DIR" EXIT

echo ""
echo "=== Running aligners ==="

echo "Running bwa-mem2..."
"$BWA_MEM2" mem -t 1 "$CHM13_INDEX" "$R1" "$R2" > "$TMP_DIR/bwa.sam" 2>/dev/null

echo "Running swiftbwa..."
"$SWIFTBWA" mem -t 1 "$CHM13_INDEX" "$R1" "$R2" > "$TMP_DIR/swift.sam" 2>/dev/null

# ─── Sort and convert to BAM ────────────────────────────────────────────────

echo ""
echo "=== Sorting and indexing BAMs ==="

for name in bwa swift; do
    echo "  $name: sort + index..."
    # Add read group if missing (GATK requires it)
    "$SAMTOOLS" addreplacerg -r '@RG\tID:bench\tSM:HG002\tPL:ILLUMINA' \
        "$TMP_DIR/${name}.sam" 2>/dev/null \
        | "$SAMTOOLS" sort -o "$TMP_DIR/${name}.bam" 2>/dev/null
    "$SAMTOOLS" index "$TMP_DIR/${name}.bam" 2>/dev/null
done

echo "BAMs ready."

# ─── GATK ValidateSamFile ───────────────────────────────────────────────────

echo ""
echo "=== GATK ValidateSamFile ==="
echo ""

for name in bwa swift; do
    tool_label=$( [[ "$name" == "bwa" ]] && echo "bwa-mem2" || echo "swiftbwa" )
    echo "--- $tool_label ---"
    "$GATK" ValidateSamFile \
        -I "$TMP_DIR/${name}.bam" \
        -R "$CHM13_REF" \
        -MODE SUMMARY \
        2>/dev/null | tee "$RESULTS_DIR/validate_${name}_${NUM_READS}.txt"
    echo ""
done

# Check for errors in swiftbwa output
SWIFT_ERRORS=$(grep -c "ERROR" "$RESULTS_DIR/validate_swift_${NUM_READS}.txt" 2>/dev/null || echo "0")
if [[ "$SWIFT_ERRORS" -gt 0 ]]; then
    echo "WARNING: swiftbwa output has $SWIFT_ERRORS GATK validation errors!"
else
    echo "swiftbwa: 0 GATK validation errors"
fi

# ─── GATK CollectAlignmentSummaryMetrics ─────────────────────────────────────

echo ""
echo "=== GATK CollectAlignmentSummaryMetrics ==="
echo ""

for name in bwa swift; do
    "$GATK" CollectAlignmentSummaryMetrics \
        -I "$TMP_DIR/${name}.bam" \
        -R "$CHM13_REF" \
        -O "$TMP_DIR/${name}_align_metrics.txt" \
        2>/dev/null
done

# Parse and compare alignment metrics
echo "Metric                      bwa-mem2         swiftbwa         diff"
echo "─────────────────────────── ──────────────── ──────────────── ────────"

compare_metric() {
    local metric="$1"
    local bwa_file="$TMP_DIR/bwa_align_metrics.txt"
    local swift_file="$TMP_DIR/swift_align_metrics.txt"

    # Find the PAIR row (or UNPAIRED for SE) and extract the metric column
    # Metrics file has a header row starting with CATEGORY, then data rows
    local bwa_val swift_val

    # Get column index for this metric
    local col_idx
    col_idx=$(head -7 "$bwa_file" | grep "^CATEGORY" | tr '\t' '\n' | grep -n "^${metric}$" | cut -d: -f1)

    if [[ -z "$col_idx" ]]; then
        return
    fi

    bwa_val=$(grep "^PAIR" "$bwa_file" | head -1 | cut -f"$col_idx")
    swift_val=$(grep "^PAIR" "$swift_file" | head -1 | cut -f"$col_idx")

    if [[ -z "$bwa_val" ]]; then
        bwa_val=$(grep "^UNPAIRED" "$bwa_file" | head -1 | cut -f"$col_idx")
        swift_val=$(grep "^UNPAIRED" "$swift_file" | head -1 | cut -f"$col_idx")
    fi

    if [[ -n "$bwa_val" && -n "$swift_val" ]]; then
        local diff
        diff=$(echo "$swift_val - $bwa_val" | bc -l 2>/dev/null || echo "N/A")
        printf "%-28s %-16s %-16s %s\n" "$metric" "$bwa_val" "$swift_val" "$diff"
    fi
}

METRICS=(
    TOTAL_READS
    PF_READS_ALIGNED
    PCT_PF_READS_ALIGNED
    PF_MISMATCH_RATE
    PF_INDEL_RATE
    MEAN_READ_LENGTH
    PCT_CHIMERAS
)

for m in "${METRICS[@]}"; do
    compare_metric "$m"
done

# Save full metrics
cp "$TMP_DIR/bwa_align_metrics.txt" "$RESULTS_DIR/align_metrics_bwa_${NUM_READS}.txt"
cp "$TMP_DIR/swift_align_metrics.txt" "$RESULTS_DIR/align_metrics_swift_${NUM_READS}.txt"

# ─── GATK CollectInsertSizeMetrics ───────────────────────────────────────────

echo ""
echo "=== GATK CollectInsertSizeMetrics ==="
echo ""

for name in bwa swift; do
    "$GATK" CollectInsertSizeMetrics \
        -I "$TMP_DIR/${name}.bam" \
        -O "$TMP_DIR/${name}_insert_metrics.txt" \
        -H "$TMP_DIR/${name}_insert_hist.pdf" \
        2>/dev/null || echo "  ($name: insert size metrics failed — may need more proper pairs)"
done

echo "Metric                      bwa-mem2         swiftbwa         diff"
echo "─────────────────────────── ──────────────── ──────────────── ────────"

compare_insert_metric() {
    local metric="$1"
    local bwa_file="$TMP_DIR/bwa_insert_metrics.txt"
    local swift_file="$TMP_DIR/swift_insert_metrics.txt"

    if [[ ! -f "$bwa_file" || ! -f "$swift_file" ]]; then
        return
    fi

    local col_idx
    col_idx=$(grep "^MEDIAN_INSERT_SIZE" "$bwa_file" | head -1 | tr '\t' '\n' | grep -n "^${metric}$" | cut -d: -f1)

    if [[ -z "$col_idx" ]]; then
        # Try header line
        col_idx=$(grep -m1 "MEDIAN_INSERT_SIZE" "$bwa_file" | tr '\t' '\n' | grep -n "^${metric}$" | cut -d: -f1)
    fi

    if [[ -z "$col_idx" ]]; then
        return
    fi

    # Data is on the line after the header
    local bwa_val swift_val
    bwa_val=$(grep -A1 "^MEDIAN_INSERT_SIZE" "$bwa_file" | tail -1 | cut -f"$col_idx")
    swift_val=$(grep -A1 "^MEDIAN_INSERT_SIZE" "$swift_file" | tail -1 | cut -f"$col_idx")

    if [[ -n "$bwa_val" && -n "$swift_val" ]]; then
        local diff
        diff=$(echo "$swift_val - $bwa_val" | bc -l 2>/dev/null || echo "N/A")
        printf "%-28s %-16s %-16s %s\n" "$metric" "$bwa_val" "$swift_val" "$diff"
    fi
}

INSERT_METRICS=(
    MEDIAN_INSERT_SIZE
    MEAN_INSERT_SIZE
    STANDARD_DEVIATION
    MEDIAN_ABSOLUTE_DEVIATION
)

for m in "${INSERT_METRICS[@]}"; do
    compare_insert_metric "$m"
done

# Save full insert metrics
if [[ -f "$TMP_DIR/bwa_insert_metrics.txt" ]]; then
    cp "$TMP_DIR/bwa_insert_metrics.txt" "$RESULTS_DIR/insert_metrics_bwa_${NUM_READS}.txt"
fi
if [[ -f "$TMP_DIR/swift_insert_metrics.txt" ]]; then
    cp "$TMP_DIR/swift_insert_metrics.txt" "$RESULTS_DIR/insert_metrics_swift_${NUM_READS}.txt"
fi

# ─── Summary ─────────────────────────────────────────────────────────────────

echo ""
echo "=== Validation Summary ==="
echo ""
echo "GATK ValidateSamFile errors:"
echo "  bwa-mem2: $(grep -c "ERROR" "$RESULTS_DIR/validate_bwa_${NUM_READS}.txt" 2>/dev/null || echo "0")"
echo "  swiftbwa: $(grep -c "ERROR" "$RESULTS_DIR/validate_swift_${NUM_READS}.txt" 2>/dev/null || echo "0")"
echo ""
echo "Results saved to: $RESULTS_DIR/"
echo "  validate_bwa_${NUM_READS}.txt"
echo "  validate_swift_${NUM_READS}.txt"
echo "  align_metrics_bwa_${NUM_READS}.txt"
echo "  align_metrics_swift_${NUM_READS}.txt"
echo "  insert_metrics_bwa_${NUM_READS}.txt"
echo "  insert_metrics_swift_${NUM_READS}.txt"
