#!/usr/bin/env python3
"""
Extract genomic coordinates of a given GFF3 feature type (or set of
gene_biotype values) into a plain BED, for use as a proximity-veto track
in scripts/score_ship_candidates.py (miRNA, lncRNA/small RNA, or a
genome-wide gene list for the 50kb/300kb radius checks - see
05_SHIP/ahmed2026_checklist_comparison.md).
"""
import argparse
import re


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--feature-types", nargs="+", required=True,
                     help="GFF3 column-3 feature types to keep, e.g. miRNA")
    ap.add_argument("--gene-biotypes", nargs="*", default=None,
                     help="optional: further filter by gene_biotype= attribute value")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    feature_types = set(args.feature_types)
    biotypes = set(args.gene_biotypes) if args.gene_biotypes else None

    n_written = 0
    with open(args.gff3, encoding="utf-8") as f, open(args.out, "w", encoding="utf-8") as out:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in feature_types:
                continue
            if biotypes is not None:
                m = re.search(r"gene_biotype=([^;]+)", fields[8])
                if not m or m.group(1) not in biotypes:
                    continue
            chrom, start, end = fields[0], int(fields[3]) - 1, int(fields[4])
            name_m = re.search(r"gene=([^;]+)", fields[8]) or re.search(r"ID=([^;]+)", fields[8])
            name = name_m.group(1) if name_m else "NA"
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")
            n_written += 1

    print(f"Wrote {n_written} features ({', '.join(feature_types)}) to {args.out}")


if __name__ == "__main__":
    main()
