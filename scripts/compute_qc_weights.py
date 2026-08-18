#!/usr/bin/env python3
"""
Aggregate per-dog QC (FRiP + TSS enrichment) into a composite weight w_i used
in the Mother Track weighted mean. Each metric is min-max normalized to [0,1]
across the 71 dogs, combined as their average, then rescaled onto
[min_weight_floor, 1.0] so no dog is ever fully zeroed out - a bad-QC dog
still contributes, just less than a good one.
"""
import argparse
import csv
import glob
import os


def minmax(values):
    lo, hi = min(values), max(values)
    if hi - lo < 1e-12:
        return [1.0 for _ in values]
    return [(v - lo) / (hi - lo) for v in values]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-dog-dir", required=True,
                    help="WORK/per_dog directory containing <sample>/qc/<sample>.qc.tsv")
    ap.add_argument("--min-weight-floor", type=float, default=0.2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    qc_files = sorted(glob.glob(os.path.join(args.per_dog_dir, "*", "qc", "*.qc.tsv")))
    if not qc_files:
        raise SystemExit(f"No QC files found under {args.per_dog_dir}")

    rows = []
    for path in qc_files:
        with open(path, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows.append(next(reader))

    samples = [r["sample"] for r in rows]
    frip = [float(r["frip"]) for r in rows]
    tss = [float(r["tss_enrichment"]) for r in rows]

    frip_norm = minmax(frip)
    tss_norm = minmax(tss)

    floor = args.min_weight_floor
    weights = []
    for fn, tn in zip(frip_norm, tss_norm):
        composite = (fn + tn) / 2.0
        w = floor + (1.0 - floor) * composite
        weights.append(w)

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("sample\tfrip\ttss_enrichment\tfrip_norm\ttss_norm\tweight\n")
        for s, fr, ts, fn, tn, w in zip(samples, frip, tss, frip_norm, tss_norm, weights):
            f.write(f"{s}\t{fr:.6f}\t{ts:.6f}\t{fn:.6f}\t{tn:.6f}\t{w:.6f}\n")

    print(f"Wrote weights for {len(samples)} dogs to {args.out}")
    print(f"weight range: {min(weights):.3f} - {max(weights):.3f}")


if __name__ == "__main__":
    main()
