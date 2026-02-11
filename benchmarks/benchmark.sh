#!/usr/bin/env bash
set -euo pipefail

# benchmark.sh — Timing benchmarks comparing bwa-mem2 and swiftbwa
# Usage: ./benchmark.sh [runs_per_config]
# Default: 3 runs per configuration

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results"

RUNS="${1:-3}"

# Tool paths
BWA_MEM2="/Users/jkane/Applications/bwa-mem2/bwa-mem2"
SWIFTBWA="$PROJECT_ROOT/.build/release/swiftbwa"
SAMTOOLS="/opt/homebrew/bin/samtools"

# Reference index
CHM13_INDEX="/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa.gz"

# Benchmark FASTQs (from existing test data)
BENCH_R1="/Users/jkane/Genomics/benchmarks/test_R1.fastq"
BENCH_R2="/Users/jkane/Genomics/benchmarks/test_R2.fastq"

# ─── Verify tools ────────────────────────────────────────────────────────────

echo "=== Verifying tools ==="
fail=0
for tool in "$BWA_MEM2" "$SWIFTBWA" "$SAMTOOLS"; do
    if [[ ! -x "$tool" ]]; then
        echo "ERROR: missing executable: $tool"
        fail=1
    fi
done
if [[ $fail -ne 0 ]]; then
    echo "Run setup.sh first to build swiftbwa."
    exit 1
fi
echo "Tools verified."

mkdir -p "$RESULTS_DIR"

# ─── Benchmark configuration ────────────────────────────────────────────────

# Test matrix: (label, reads, mode, r1, r2, thread_counts)
# Thread counts tested for each config
THREAD_COUNTS=(1 4 8)

# CSV header
CSV="$RESULTS_DIR/benchmark.csv"
echo "tool,dataset,reads,mode,threads,run,wall_sec,peak_rss_kb,mapped_pct,reads_per_sec" > "$CSV"

# ─── Timing function ────────────────────────────────────────────────────────

run_benchmark() {
    local tool_name="$1"
    local tool_cmd="$2"
    local dataset="$3"
    local num_reads="$4"
    local mode="$5"
    local threads="$6"
    local run_num="$7"
    local r1="$8"
    local r2="${9:-}"

    local tmp_sam
    tmp_sam=$(mktemp /tmp/bench_XXXXXX.sam)
    local tmp_time
    tmp_time=$(mktemp /tmp/bench_time_XXXXXX.txt)

    local cmd
    if [[ "$mode" == "SE" ]]; then
        cmd="$tool_cmd mem -t $threads $CHM13_INDEX $r1"
    else
        cmd="$tool_cmd mem -t $threads $CHM13_INDEX $r1 $r2"
    fi

    # Use /usr/bin/time -l for wall time + peak RSS (macOS)
    /usr/bin/time -l bash -c "$cmd > $tmp_sam 2>/dev/null" 2>"$tmp_time" || true

    # Parse wall time (macOS format: "N.NN real ..." or GNU: "... elapsed")
    local wall_sec
    wall_sec=$(awk '/real/{print $1; exit}' "$tmp_time" 2>/dev/null || echo "0")
    if [[ -z "$wall_sec" || "$wall_sec" == "0" ]]; then
        # Try macOS format: first line "    N.NN real    ..."
        wall_sec=$(head -1 "$tmp_time" | awk '{print $1}' 2>/dev/null || echo "0")
    fi

    # Parse peak RSS (macOS: "maximum resident set size" in bytes)
    local peak_rss_kb
    peak_rss_kb=$(awk '/maximum resident set size/{print int($1/1024); exit}' "$tmp_time" 2>/dev/null || echo "0")
    if [[ "$peak_rss_kb" == "0" ]]; then
        # Alternative: some macOS versions report differently
        peak_rss_kb=$(grep -i "maximum resident" "$tmp_time" | awk '{print int($1/1024)}' 2>/dev/null || echo "0")
    fi

    # Get mapped percentage
    local mapped_pct="0"
    if [[ -s "$tmp_sam" ]]; then
        local flagstat
        flagstat=$("$SAMTOOLS" flagstat "$tmp_sam" 2>/dev/null || echo "")
        mapped_pct=$(echo "$flagstat" | awk '/mapped \(/{gsub(/[()%]/,""); print $5; exit}' 2>/dev/null || echo "0")
    fi

    # Calculate reads/sec
    local reads_per_sec
    if (( $(echo "$wall_sec > 0" | bc -l 2>/dev/null || echo 0) )); then
        reads_per_sec=$(echo "scale=1; $num_reads / $wall_sec" | bc -l 2>/dev/null || echo "0")
    else
        reads_per_sec="0"
    fi

    echo "$tool_name,$dataset,$num_reads,$mode,$threads,$run_num,$wall_sec,$peak_rss_kb,$mapped_pct,$reads_per_sec" >> "$CSV"

    printf "    %-10s t=%-2s run=%d  %7.2fs  %6dMB  mapped=%s%%  %s reads/s\n" \
        "$tool_name" "$threads" "$run_num" "$wall_sec" "$((peak_rss_kb / 1024))" "$mapped_pct" "$reads_per_sec"

    rm -f "$tmp_sam" "$tmp_time"
}

# ─── Run benchmarks ─────────────────────────────────────────────────────────

echo ""
echo "=== Benchmark: bwa-mem2 vs swiftbwa ==="
echo "Runs per config: $RUNS (run 1 = warmup, median of runs 2+ reported)"
echo ""

# Dataset definitions: label num_reads mode r1 r2
declare -a DATASETS
DATASETS=(
    "10K_PE|10000|PE|$DATA_DIR/hg002_10000_R1.fq|$DATA_DIR/hg002_10000_R2.fq"
    "100K_PE|100000|PE|$DATA_DIR/hg002_100000_R1.fq|$DATA_DIR/hg002_100000_R2.fq"
    "100K_SE|100000|SE|$DATA_DIR/hg002_100000_R1.fq|"
    "500K_PE|500000|PE|$BENCH_R1|$BENCH_R2"
)

for ds_entry in "${DATASETS[@]}"; do
    IFS='|' read -r label num_reads mode r1 r2 <<< "$ds_entry"

    # Check if input files exist
    if [[ ! -f "$r1" ]]; then
        echo "Skipping $label: $r1 not found (run setup.sh first)"
        continue
    fi
    if [[ "$mode" == "PE" && -n "$r2" && ! -f "$r2" ]]; then
        echo "Skipping $label: $r2 not found (run setup.sh first)"
        continue
    fi

    echo "--- $label ($num_reads reads, $mode) ---"

    for T in "${THREAD_COUNTS[@]}"; do
        for run in $(seq 1 "$RUNS"); do
            run_benchmark "bwa-mem2" "$BWA_MEM2" "$label" "$num_reads" "$mode" "$T" "$run" "$r1" "$r2"
            run_benchmark "swiftbwa" "$SWIFTBWA" "$label" "$num_reads" "$mode" "$T" "$run" "$r1" "$r2"
        done
    done
    echo ""
done

# ─── Summary table ───────────────────────────────────────────────────────────

echo "=== Summary (median of non-warmup runs) ==="
echo ""
printf "%-12s %-10s %-4s %3s  %10s %10s  %8s %8s  %7s\n" \
    "Dataset" "Tool" "Mode" "T" "Wall(s)" "Reads/s" "Mapped%" "RSS(MB)" "Speedup"
printf "%s\n" "$(printf '%.0s-' {1..90})"

# Process CSV to compute median of runs 2+
awk -F',' '
NR == 1 { next }
{
    key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4 SUBSEP $5
    if ($6 > 1) {
        walls[key][++n[key]] = $7
        rps[key][n[key]] = $10
        mapped[key] = $9
        rss[key] = $8
    }
}
END {
    for (key in n) {
        split(key, parts, SUBSEP)
        tool = parts[1]; ds = parts[2]; reads = parts[3]; mode = parts[4]; threads = parts[5]
        cnt = n[key]
        # Sort walls for median
        for (i = 1; i <= cnt; i++)
            for (j = i+1; j <= cnt; j++)
                if (walls[key][i] > walls[key][j]) {
                    tmp = walls[key][i]; walls[key][i] = walls[key][j]; walls[key][j] = tmp
                    tmp = rps[key][i]; rps[key][i] = rps[key][j]; rps[key][j] = tmp
                }
        mid = int((cnt+1)/2)
        printf "%s,%s,%s,%s,%s,%.2f,%.1f,%s,%s\n", \
            ds, tool, mode, threads, reads, walls[key][mid], rps[key][mid], mapped[key], rss[key]
    }
}' "$CSV" | sort -t',' -k1,1 -k4,4n -k2,2 | while IFS=',' read -r ds tool mode threads reads wall rps mapped_pct rss; do
    # Look up bwa-mem2 wall time for speedup calculation
    printf "%-12s %-10s %-4s %3s  %10.2f %10.1f  %7s%% %7dMB\n" \
        "$ds" "$tool" "$mode" "$threads" "$wall" "$rps" "$mapped_pct" "$((rss / 1024))"
done

echo ""
echo "Full results: $CSV"
echo "Total rows: $(( $(wc -l < "$CSV") - 1 ))"
