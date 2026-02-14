#!/usr/bin/env python3
"""Supplementary alignment analysis: compare swiftbwa vs bwa-mem2.

Usage: python3 supp_analysis.py --swift-sam PATH --bwa-sam PATH [--limit N]
"""

import argparse
import sys
from collections import defaultdict


def parse_sam(sam_path):
    """Parse SAM into dict: readname -> list of records (primaries + supps)."""
    reads = defaultdict(list)
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            flag = int(fields[1])
            if flag & 0x100:
                continue  # skip secondary
            key = fields[0]
            if flag & 0x40:
                key += "/1"
            elif flag & 0x80:
                key += "/2"
            reads[key].append({
                "flag": flag,
                "rname": fields[2],
                "pos": int(fields[3]),
                "mapq": int(fields[4]),
                "cigar": fields[5],
                "is_supp": bool(flag & 0x800),
                "is_unmapped": bool(flag & 0x4),
                "is_reverse": bool(flag & 0x10),
            })
    return reads


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--swift-sam", required=True)
    parser.add_argument("--bwa-sam", required=True)
    parser.add_argument("--limit", type=int, default=20,
                        help="Max examples to print")
    args = parser.parse_args()

    print("Parsing bwa-mem2 SAM...")
    bwa = parse_sam(args.bwa_sam)
    print("Parsing swiftbwa SAM...")
    swift = parse_sam(args.swift_sam)

    # Count supplementaries per read
    bwa_supp_counts = defaultdict(int)
    swift_supp_counts = defaultdict(int)
    for key, recs in bwa.items():
        bwa_supp_counts[key] = sum(1 for r in recs if r["is_supp"])
    for key, recs in swift.items():
        swift_supp_counts[key] = sum(1 for r in recs if r["is_supp"])

    # Categorize
    swift_only_supp = []  # reads with supps in swift but not bwa
    bwa_only_supp = []    # reads with supps in bwa but not swift
    both_supp = []        # reads with supps in both
    swift_more_supp = []  # reads with more supps in swift

    all_keys = set(bwa_supp_counts) | set(swift_supp_counts)
    for key in sorted(all_keys):
        sc = swift_supp_counts.get(key, 0)
        bc = bwa_supp_counts.get(key, 0)
        if sc > 0 and bc == 0:
            swift_only_supp.append((key, sc))
        elif bc > 0 and sc == 0:
            bwa_only_supp.append((key, bc))
        elif sc > 0 and bc > 0:
            both_supp.append((key, sc, bc))
            if sc > bc:
                swift_more_supp.append((key, sc, bc))

    print()
    print("=== Supplementary Distribution ===")
    print(f"  Reads with supps in swift only:  {len(swift_only_supp)}")
    print(f"  Reads with supps in bwa only:    {len(bwa_only_supp)}")
    print(f"  Reads with supps in both:        {len(both_supp)}")
    print(f"  Reads with more supps in swift:  {len(swift_more_supp)}")
    print()

    # Analyze swift-only supp reads
    print("=== Swift-Only Supp Reads: Characteristics ===")
    mapq_dist = defaultdict(int)
    supp_mapq_dist = defaultdict(int)
    primary_positions_match = 0
    primary_positions_differ = 0

    for key, sc in swift_only_supp:
        sw_recs = swift[key]
        bw_recs = bwa.get(key, [])
        # Look at the primary
        sw_pri = [r for r in sw_recs if not r["is_supp"]]
        sw_supps = [r for r in sw_recs if r["is_supp"]]
        bw_pri = [r for r in bw_recs if not r["is_supp"]]

        for sp in sw_supps:
            supp_mapq_dist[sp["mapq"]] += 1

        if sw_pri and bw_pri:
            sp = sw_pri[0]
            bp = bw_pri[0]
            mapq_dist[sp["mapq"]] += 1
            if sp["rname"] == bp["rname"] and sp["pos"] == bp["pos"]:
                primary_positions_match += 1
            else:
                primary_positions_differ += 1

    print(f"  Primary position matches bwa:    {primary_positions_match}")
    print(f"  Primary position differs:        {primary_positions_differ}")
    print()

    # Supp MAPQ distribution
    print("  Supplementary MAPQ distribution:")
    for mq in sorted(supp_mapq_dist.keys()):
        print(f"    MAPQ {mq:3d}: {supp_mapq_dist[mq]}")
    print()

    # Primary MAPQ distribution for reads with swift-only supps
    print("  Primary MAPQ distribution (swift):")
    for mq in sorted(mapq_dist.keys()):
        print(f"    MAPQ {mq:3d}: {mapq_dist[mq]}")
    print()

    # Show examples
    print(f"=== Examples: Swift has supps, BWA does not (first {args.limit}) ===")
    count = 0
    for key, sc in swift_only_supp:
        if count >= args.limit:
            break
        sw_recs = swift[key]
        bw_recs = bwa.get(key, [])
        print(f"\n  Read: {key}  (swift supps: {sc})")
        print("    SWIFT:")
        for r in sw_recs:
            tag = "SUPP" if r["is_supp"] else "PRI "
            strand = "-" if r["is_reverse"] else "+"
            print(f"      {tag} {r['rname']}:{r['pos']} {strand} MAPQ={r['mapq']} CIGAR={r['cigar'][:60]}")
        print("    BWA:")
        for r in bw_recs:
            tag = "SUPP" if r["is_supp"] else "PRI "
            strand = "-" if r["is_reverse"] else "+"
            print(f"      {tag} {r['rname']}:{r['pos']} {strand} MAPQ={r['mapq']} CIGAR={r['cigar'][:60]}")
        count += 1


if __name__ == "__main__":
    main()
