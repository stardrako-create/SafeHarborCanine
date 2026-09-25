#!/usr/bin/env python3
"""
One-time v1.20.0 backfill: add the `lambda` column to the 76 dogs'
already-existing qc.tsv files, without rerunning the per-sample Snakemake
pipeline (that pipeline is unchanged in structure/results - lambda is
computed from the same raw.bw values it would have used anyway, just
computed now instead of at the original run time). Reads the .raw.bw
header (fast, no full-file scan) for each dog and rewrites its qc.tsv with
lambda appended.

Run once against the unified 76-dog sample list (71 Jin2024 + 5 Ehsan).
After this, build_mother_track_v2.py never needs to open raw.bw again.
"""
import argparse
import csv
import os

import pyBigWig


def compute_lambda(raw_bw_path):
    bw = pyBigWig.open(raw_bw_path)
    header = bw.header()
    genome_bp = sum(bw.chroms().values())
    bw.close()
    return max(header.get("sumData", 0.0) / genome_bp, 1e-6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, nargs="+")
    ap.add_argument("--per-dog-dir", required=True)
    args = ap.parse_args()

    for s in args.samples:
        qc_path = os.path.join(args.per_dog_dir, s, "qc", f"{s}.qc.tsv")
        raw_bw_path = os.path.join(args.per_dog_dir, s, "bigwig", f"{s}.raw.bw")

        with open(qc_path, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = reader.fieldnames
            row = next(reader)

        if "lambda" in fieldnames:
            print(f"{s}: qc.tsv already has lambda, skipping")
            continue

        lam = compute_lambda(raw_bw_path)
        row["lambda"] = f"{lam:.6f}"
        new_fieldnames = fieldnames + ["lambda"]

        with open(qc_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=new_fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerow(row)

        print(f"{s}: lambda={lam:.4f}")


if __name__ == "__main__":
    main()
