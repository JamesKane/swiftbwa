#!/usr/bin/env python3
"""MAPQ overconfidence analysis between bwa-mem2 and swiftbwa.

Usage: python3 mapq_analysis.py [--reads N] [--threads N] [--swift-sam PATH] [--bwa-sam PATH]

By default, runs swiftbwa on hg002_{reads} data and uses pre-computed bwa-mem2 SAM.
If --swift-sam is given, skips running swiftbwa.
If --bwa-sam is given, uses that instead of the pre-computed reference.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(SCRIPT_DIR, "data")

SWIFTBWA = os.path.join(PROJECT_ROOT, ".build", "release", "swiftbwa")
BWA_MEM2 = "/Users/jkane/Applications/bwa-mem2/bwa-mem2"
CHM13_INDEX = "/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa.gz"


def parse_sam_primaries(sam_path):
    """Parse primary alignments from SAM, keyed by readname/1 or /2."""
    records = {}
    supp_count = 0
    mapped_count = 0
    proper_count = 0
    total = 0
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            if flag & 0x800:
                supp_count += 1
                continue
            if flag & 0x100:
                continue
            # Primary alignment
            total += 1
            if not (flag & 0x4):
                mapped_count += 1
            if flag & 0x2:
                proper_count += 1

            key = fields[0]
            if flag & 0x40:
                key += "/1"
            elif flag & 0x80:
                key += "/2"

            records[key] = {
                "flag": flag,
                "rname": fields[2],
                "pos": int(fields[3]),
                "mapq": int(fields[4]),
                "cigar": fields[5],
            }
    return records, total, mapped_count, proper_count, supp_count


def run_analysis(bwa_sam, swift_sam):
    print("Parsing bwa-mem2 SAM...")
    bwa, bwa_total, bwa_mapped, bwa_proper, bwa_supp = parse_sam_primaries(bwa_sam)
    print("Parsing swiftbwa SAM...")
    swift, sw_total, sw_mapped, sw_proper, sw_supp = parse_sam_primaries(swift_sam)

    both_mapped = 0
    pos_exact = 0
    pos_mismatch = 0
    overconfident = 0
    oc_mapq_sum = 0
    both_mapq0 = 0
    high_conf_mismatch = 0
    mismatch_sw_mapq_sum = 0
    mismatch_bwa_mapq_sum = 0
    mismatch_count = 0

    # MAPQ distribution for mismatched reads
    mapq_exact = 0
    mapq_total = 0
    mapq_diff_sum = 0

    for key, sr in swift.items():
        br = bwa.get(key)
        if br is None:
            continue
        if (br["flag"] & 0x4) or (sr["flag"] & 0x4):
            continue

        both_mapped += 1
        mapq_total += 1
        mapq_diff_sum += abs(br["mapq"] - sr["mapq"])
        if br["mapq"] == sr["mapq"]:
            mapq_exact += 1

        same_pos = br["rname"] == sr["rname"] and br["pos"] == sr["pos"]
        if same_pos:
            pos_exact += 1
        else:
            pos_mismatch += 1
            mismatch_count += 1
            sm = sr["mapq"]
            bm = br["mapq"]
            mismatch_sw_mapq_sum += sm
            mismatch_bwa_mapq_sum += bm

            if sm >= 50 and bm <= 5:
                overconfident += 1
                oc_mapq_sum += sm
            if sm == 0 and bm == 0:
                both_mapq0 += 1
            if sm >= 30 and bm >= 30:
                high_conf_mismatch += 1

    print()
    print("=== Results ===")
    print()
    print("-- Mapping Rate --")
    print(f"  bwa-mem2:  {bwa_mapped} / {bwa_total} ({100*bwa_mapped/bwa_total:.1f}%)")
    print(f"  swiftbwa:  {sw_mapped} / {sw_total} ({100*sw_mapped/sw_total:.1f}%)")
    print()
    print("-- Position Concordance --")
    print(f"  Both mapped:     {both_mapped}")
    print(f"  Exact match:     {pos_exact} ({100*pos_exact/both_mapped:.2f}%)")
    print(f"  Mismatch:        {pos_mismatch} ({100*pos_mismatch/both_mapped:.2f}%)")
    print()
    print("-- MAPQ Agreement (both mapped) --")
    print(f"  Exact match:     {mapq_exact} / {mapq_total} ({100*mapq_exact/mapq_total:.1f}%)")
    print(f"  Mean abs diff:   {mapq_diff_sum/mapq_total:.2f}")
    print()
    print("-- Overconfident MAPQ (swift>=50, bwa<=5, diff pos) --")
    print(f"  Count:           {overconfident}")
    if overconfident > 0:
        print(f"  Mean swift MAPQ: {oc_mapq_sum/overconfident:.1f}")
    print()
    print("-- Mismatched Position Reads --")
    if mismatch_count > 0:
        print(f"  Count:           {mismatch_count}")
        print(f"  Mean swift MAPQ: {mismatch_sw_mapq_sum/mismatch_count:.1f}")
        print(f"  Mean bwa MAPQ:   {mismatch_bwa_mapq_sum/mismatch_count:.1f}")
    print(f"  Both MAPQ=0:     {both_mapq0}")
    print(f"  Both MAPQ>=30:   {high_conf_mismatch}")
    print()
    print("-- Proper Pair Rate --")
    print(f"  bwa-mem2:  {100*bwa_proper/bwa_total:.1f}%")
    print(f"  swiftbwa:  {100*sw_proper/sw_total:.1f}%")
    print()
    print("-- Supplementary Alignments --")
    print(f"  bwa-mem2:  {bwa_supp}")
    print(f"  swiftbwa:  {sw_supp}")
    print()


def main():
    parser = argparse.ArgumentParser(description="MAPQ overconfidence analysis")
    parser.add_argument("--reads", type=int, default=100000)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--swift-sam", help="Pre-computed swiftbwa SAM (skip running)")
    parser.add_argument("--bwa-sam", help="Pre-computed bwa-mem2 SAM")
    parser.add_argument("--r1", help="Read 1 FASTQ path")
    parser.add_argument("--r2", help="Read 2 FASTQ path")
    args = parser.parse_args()

    r1 = args.r1 or os.path.join(DATA_DIR, f"hg002_{args.reads}_R1.fq")
    r2 = args.r2 or os.path.join(DATA_DIR, f"hg002_{args.reads}_R2.fq")
    bwa_sam = args.bwa_sam or os.path.join(DATA_DIR, f"bwamem2_pe_{args.reads}.sam")

    print(f"=== MAPQ Overconfidence Analysis ===")
    print(f"Reads: {args.reads}  Threads: {args.threads}")
    print()

    # Run or reuse swiftbwa
    if args.swift_sam:
        swift_sam = args.swift_sam
    else:
        if not os.path.isfile(SWIFTBWA):
            print(f"ERROR: build swiftbwa first: {SWIFTBWA}")
            sys.exit(1)
        for f in [r1, r2]:
            if not os.path.isfile(f):
                print(f"ERROR: missing: {f}")
                sys.exit(1)

        swift_sam = tempfile.mktemp(suffix=".sam", prefix="mapq_swift_")
        print("Running swiftbwa...")
        with open(swift_sam, "w") as out:
            subprocess.run(
                [SWIFTBWA, "mem", "-t", str(args.threads), CHM13_INDEX, r1, r2],
                stdout=out, stderr=subprocess.DEVNULL, check=True,
            )
        print("Done.")
        print()

    if not os.path.isfile(bwa_sam):
        print(f"ERROR: missing bwa-mem2 SAM: {bwa_sam}")
        sys.exit(1)

    try:
        run_analysis(bwa_sam, swift_sam)
    finally:
        if not args.swift_sam and os.path.exists(swift_sam):
            os.unlink(swift_sam)


if __name__ == "__main__":
    main()
