#!/usr/bin/env python3
"""Deep supplementary analysis: check SA/XA tags and query overlap patterns."""

import argparse
import sys
from collections import defaultdict


def parse_sam_full(sam_path):
    """Parse SAM into dict: readname -> list of all records with tags."""
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

            # Parse optional tags
            tags = {}
            for t in fields[11:]:
                parts = t.split(":", 2)
                if len(parts) >= 3:
                    tags[parts[0]] = parts[2]

            reads[key].append({
                "flag": flag,
                "rname": fields[2],
                "pos": int(fields[3]),
                "mapq": int(fields[4]),
                "cigar": fields[5],
                "is_supp": bool(flag & 0x800),
                "is_unmapped": bool(flag & 0x4),
                "is_reverse": bool(flag & 0x10),
                "SA": tags.get("SA", ""),
                "XA": tags.get("XA", ""),
            })
    return reads


def count_xa_entries(xa_tag):
    """Count number of XA tag entries."""
    if not xa_tag:
        return 0
    return xa_tag.count(";")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--swift-sam", required=True)
    parser.add_argument("--bwa-sam", required=True)
    parser.add_argument("--limit", type=int, default=30)
    args = parser.parse_args()

    print("Parsing...")
    bwa = parse_sam_full(args.bwa_sam)
    swift = parse_sam_full(args.swift_sam)

    # Find reads where swift has supplementaries but bwa doesn't
    swift_only_supp = []
    for key in swift:
        sw_supps = sum(1 for r in swift[key] if r["is_supp"])
        bw_supps = sum(1 for r in bwa.get(key, []) if r["is_supp"]) if key in bwa else 0
        if sw_supps > 0 and bw_supps == 0:
            swift_only_supp.append(key)

    print(f"Swift-only supplementary reads: {len(swift_only_supp)}")

    # Categorize: does bwa have XA tag (reporting these as secondary)?
    bwa_has_xa = 0
    bwa_no_xa = 0
    bwa_xa_count_dist = defaultdict(int)
    swift_sa_count_dist = defaultdict(int)
    swift_xa_count_dist = defaultdict(int)

    # Check if swift supp position appears in bwa's XA tag
    supp_in_bwa_xa = 0
    supp_not_in_bwa_xa = 0

    for key in swift_only_supp:
        sw_recs = swift[key]
        bw_recs = bwa.get(key, [])

        # Swift supplementary info
        sw_supps = [r for r in sw_recs if r["is_supp"]]
        sw_pri = [r for r in sw_recs if not r["is_supp"]]

        sa_count = sum(1 for r in sw_pri if r.get("SA"))
        swift_sa_count_dist[sa_count] += 1
        for r in sw_pri:
            swift_xa_count_dist[count_xa_entries(r["XA"])] += 1

        # BWA primary XA tag
        bw_pri = [r for r in bw_recs if not r["is_supp"]]
        if bw_pri:
            xa = bw_pri[0].get("XA", "")
            xa_count = count_xa_entries(xa)
            bwa_xa_count_dist[xa_count] += 1
            if xa_count > 0:
                bwa_has_xa += 1
            else:
                bwa_no_xa += 1

            # Check if any swift supplementary position appears in bwa XA
            for supp in sw_supps:
                supp_chr = supp["rname"]
                supp_pos = supp["pos"]
                supp_strand = "-" if supp["is_reverse"] else "+"
                # XA format: chr,+pos,CIGAR,NM;
                search = f"{supp_chr},{supp_strand}{supp_pos},"
                if search in xa:
                    supp_in_bwa_xa += 1
                else:
                    supp_not_in_bwa_xa += 1

    print()
    print("=== BWA XA tag for reads with swift-only supps ===")
    print(f"  BWA primary has XA:    {bwa_has_xa}")
    print(f"  BWA primary no XA:     {bwa_no_xa}")
    print(f"  Swift supp pos in bwa XA:  {supp_in_bwa_xa}")
    print(f"  Swift supp pos NOT in bwa XA: {supp_not_in_bwa_xa}")
    print()

    print("BWA XA entry count distribution:")
    for c in sorted(bwa_xa_count_dist.keys()):
        print(f"    {c} entries: {bwa_xa_count_dist[c]}")
    print()

    print("Swift SA tag distribution (primary):")
    for c in sorted(swift_sa_count_dist.keys()):
        print(f"    {c} SA tags: {swift_sa_count_dist[c]}")
    print()

    # Show detailed examples
    print(f"=== Detailed examples (first {args.limit}) ===")
    count = 0
    for key in swift_only_supp:
        if count >= args.limit:
            break
        sw_recs = swift[key]
        bw_recs = bwa.get(key, [])

        print(f"\n  Read: {key}")
        print("    SWIFT:")
        for r in sw_recs:
            tag = "SUPP" if r["is_supp"] else "PRI "
            strand = "-" if r["is_reverse"] else "+"
            xa_n = count_xa_entries(r["XA"])
            sa = f" SA={r['SA'][:80]}" if r["SA"] else ""
            xa_info = f" XA({xa_n})" if xa_n else ""
            print(f"      {tag} {r['rname']}:{r['pos']} {strand} MQ={r['mapq']} {r['cigar'][:50]}{sa}{xa_info}")
        print("    BWA:")
        for r in bw_recs:
            tag = "SUPP" if r["is_supp"] else "PRI "
            strand = "-" if r["is_reverse"] else "+"
            xa_n = count_xa_entries(r["XA"])
            sa = f" SA={r['SA'][:80]}" if r["SA"] else ""
            xa_info = f" XA({xa_n})" if xa_n else ""
            print(f"      {tag} {r['rname']}:{r['pos']} {strand} MQ={r['mapq']} {r['cigar'][:50]}{sa}{xa_info}")
        count += 1


if __name__ == "__main__":
    main()
