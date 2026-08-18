#!/usr/bin/env python3
"""Extract one TSS coordinate per gene feature from a GFF3 into a BED file."""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    n = 0
    with open(args.gff, encoding="utf-8") as fin, open(args.out, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            chrom, start, end, strand = cols[0], int(cols[3]), int(cols[4]), cols[6]
            name = "gene"
            for field in cols[8].split(";"):
                if field.startswith("Name="):
                    name = field[len("Name="):]
                    break
            if strand == "-":
                tss = end - 1  # 0-based BED start
            else:
                tss = start - 1
            if tss < 0:
                continue
            fout.write(f"{chrom}\t{tss}\t{tss + 1}\t{name}\t.\t{strand}\n")
            n += 1

    print(f"Wrote {n} TSS records to {args.out}")


if __name__ == "__main__":
    main()
