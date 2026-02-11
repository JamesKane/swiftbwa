#!/usr/bin/env bash
set -euo pipefail

# setup.sh — Subsample HG002 reads and generate bwa-mem2 gold standards
# for correctness testing against chm13v2.0 reference.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
TEST_DATA_DIR="$PROJECT_ROOT/TestData"

# Tool paths
BWA_MEM2="/Users/jkane/Applications/bwa-mem2/bwa-mem2"
SAMTOOLS="/opt/homebrew/bin/samtools"

# Reference and index
CHM13_INDEX="/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa.gz"

# HG002 input FASTQs
HG002_R1="/Users/jkane/Genomics/HG002/2A1_CGATGT_L001_R1_001.fastq.gz"
HG002_R2="/Users/jkane/Genomics/HG002/2A1_CGATGT_L001_R2_001.fastq.gz"

# Subsample sizes
SIZES=(100 1000 10000 100000)

# ─── Verify tools ───────────────────────────────────────────────────────────

echo "=== Verifying tools and data ==="

fail=0
for tool in "$BWA_MEM2" "$SAMTOOLS"; do
    if [[ ! -x "$tool" ]]; then
        echo "ERROR: missing executable: $tool"
        fail=1
    fi
done

for f in "$CHM13_INDEX.bwt.2bit.64" "$CHM13_INDEX.pac" "$CHM13_INDEX.ann" "$CHM13_INDEX.amb"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing index file: $f"
        fail=1
    fi
done

for f in "$HG002_R1" "$HG002_R2"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: missing FASTQ: $f"
        fail=1
    fi
done

if [[ $fail -ne 0 ]]; then
    echo "Setup failed — fix missing dependencies above."
    exit 1
fi

echo "All tools and data verified."

# ─── Build swiftbwa ─────────────────────────────────────────────────────────

echo ""
echo "=== Building swiftbwa (release) ==="
cd "$PROJECT_ROOT"
swift build -c release 2>&1 | tail -5
SWIFTBWA="$PROJECT_ROOT/.build/release/swiftbwa"
if [[ ! -x "$SWIFTBWA" ]]; then
    echo "ERROR: swiftbwa binary not found after build"
    exit 1
fi
echo "Built: $SWIFTBWA"

# ─── Create output directories ──────────────────────────────────────────────

mkdir -p "$DATA_DIR"

# ─── Subsample reads ────────────────────────────────────────────────────────

echo ""
echo "=== Subsampling HG002 reads ==="

subsample() {
    local n=$1
    local r1_out="$2"
    local r2_out="$3"
    local lines=$((n * 4))

    if [[ -f "$r1_out" && -f "$r2_out" ]]; then
        local existing_lines
        existing_lines=$(wc -l < "$r1_out" | tr -d ' ')
        if [[ "$existing_lines" -eq "$lines" ]]; then
            echo "  $n reads: already exists, skipping"
            return
        fi
    fi

    echo "  $n reads: subsampling..."
    # head closes pipe early → gunzip gets SIGPIPE (exit 141); that's expected
    (gunzip -c "$HG002_R1" || true) | head -n "$lines" > "$r1_out"
    (gunzip -c "$HG002_R2" || true) | head -n "$lines" > "$r2_out"

    # Verify line counts
    local r1_lines r2_lines
    r1_lines=$(wc -l < "$r1_out" | tr -d ' ')
    r2_lines=$(wc -l < "$r2_out" | tr -d ' ')
    if [[ "$r1_lines" -ne "$lines" || "$r2_lines" -ne "$lines" ]]; then
        echo "  WARNING: expected $lines lines, got R1=$r1_lines R2=$r2_lines"
    fi
}

# 100 reads → TestData/ (committed, for Swift tests)
subsample 100 "$TEST_DATA_DIR/hg002_chm13_pe_100_R1.fq" "$TEST_DATA_DIR/hg002_chm13_pe_100_R2.fq"
# Also copy R1 as SE input
cp "$TEST_DATA_DIR/hg002_chm13_pe_100_R1.fq" "$TEST_DATA_DIR/hg002_chm13_se_100.fq"

# Larger sizes → benchmarks/data/ (not committed)
for n in "${SIZES[@]}"; do
    if [[ "$n" -eq 100 ]]; then
        # Already in TestData; also copy to data/ for consistency
        cp "$TEST_DATA_DIR/hg002_chm13_pe_100_R1.fq" "$DATA_DIR/hg002_100_R1.fq"
        cp "$TEST_DATA_DIR/hg002_chm13_pe_100_R2.fq" "$DATA_DIR/hg002_100_R2.fq"
    else
        subsample "$n" "$DATA_DIR/hg002_${n}_R1.fq" "$DATA_DIR/hg002_${n}_R2.fq"
    fi
done

echo "Subsampling complete."

# ─── Generate bwa-mem2 gold standards ────────────────────────────────────────

echo ""
echo "=== Generating bwa-mem2 gold standards ==="

# 100-read gold standards for Swift tests (headerless, committed to TestData/)
echo "  SE 100 reads (TestData)..."
"$BWA_MEM2" mem -t 1 "$CHM13_INDEX" \
    "$TEST_DATA_DIR/hg002_chm13_se_100.fq" 2>/dev/null \
    | grep -v '^@' > "$TEST_DATA_DIR/hg002_chm13_se_100.sam"

echo "  PE 100 reads (TestData)..."
"$BWA_MEM2" mem -t 1 "$CHM13_INDEX" \
    "$TEST_DATA_DIR/hg002_chm13_pe_100_R1.fq" \
    "$TEST_DATA_DIR/hg002_chm13_pe_100_R2.fq" 2>/dev/null \
    | grep -v '^@' > "$TEST_DATA_DIR/hg002_chm13_pe_100.sam"

# Larger gold standards in benchmarks/data/ (not committed)
for n in "${SIZES[@]}"; do
    local_r1="$DATA_DIR/hg002_${n}_R1.fq"
    local_r2="$DATA_DIR/hg002_${n}_R2.fq"

    if [[ ! -f "$local_r1" ]]; then
        echo "  Skipping $n (no R1 file)"
        continue
    fi

    echo "  PE $n reads (benchmarks/data)..."
    "$BWA_MEM2" mem -t 1 "$CHM13_INDEX" "$local_r1" "$local_r2" 2>/dev/null \
        | grep -v '^@' > "$DATA_DIR/bwamem2_pe_${n}.sam"
done

echo "Gold standards generated."

# ─── Summary ─────────────────────────────────────────────────────────────────

echo ""
echo "=== Summary ==="
echo ""
echo "TestData/ (committed, for Swift tests):"
for f in "$TEST_DATA_DIR"/hg002_chm13_*; do
    if [[ -f "$f" ]]; then
        local_lines=$(wc -l < "$f" | tr -d ' ')
        local_size=$(du -h "$f" | cut -f1)
        echo "  $(basename "$f"): $local_lines lines, $local_size"
    fi
done

echo ""
echo "benchmarks/data/ (not committed):"
for f in "$DATA_DIR"/*; do
    if [[ -f "$f" ]]; then
        local_lines=$(wc -l < "$f" | tr -d ' ')
        local_size=$(du -h "$f" | cut -f1)
        echo "  $(basename "$f"): $local_lines lines, $local_size"
    fi
done

# Quick sanity check with samtools
echo ""
echo "=== Sanity check (samtools flagstat on 100-read PE gold) ==="
"$SAMTOOLS" flagstat <("$BWA_MEM2" mem -t 1 "$CHM13_INDEX" \
    "$TEST_DATA_DIR/hg002_chm13_pe_100_R1.fq" \
    "$TEST_DATA_DIR/hg002_chm13_pe_100_R2.fq" 2>/dev/null)

echo ""
echo "Setup complete!"
