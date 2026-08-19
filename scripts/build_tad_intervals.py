#!/usr/bin/env python3
"""
Derive TAD (topologically associating domain) intervals from the boundary
calls in 04_tracks_processadas/.../HiC/tad_boundaries.bed: adjacent/
overlapping 25kb boundary bins are merged into single boundary blocks per
chromosome, then a TAD is the span between the trailing edge of one merged
boundary block and the leading edge of the next. Used to check whether a
candidate's whole TAD (not just its two flanking genes) contains a cancer-
related gene - criterion 8 in Ahmed et al. 2026 - see
05_SHIP/ahmed2026_checklist_comparison.md.
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--boundaries-bed", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    by_chrom = {}
    with open(args.boundaries_bed, encoding="utf-8") as f:
        for line in f:
            chrom, start, end = line.rstrip("\n").split("\t")[:3]
            by_chrom.setdefault(chrom, []).append((int(start), int(end)))

    n_tads = 0
    with open(args.out, "w", encoding="utf-8") as out:
        for chrom, intervals in by_chrom.items():
            intervals.sort()
            merged = []
            for s, e in intervals:
                if merged and s <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((s, e))
            for i in range(len(merged) - 1):
                tad_start, tad_end = merged[i][1], merged[i + 1][0]
                if tad_end > tad_start:
                    out.write(f"{chrom}\t{tad_start}\t{tad_end}\n")
                    n_tads += 1

    print(f"Wrote {n_tads} TAD intervals to {args.out}")


if __name__ == "__main__":
    main()
