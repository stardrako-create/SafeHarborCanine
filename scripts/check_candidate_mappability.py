#!/usr/bin/env python3
"""
Self-mappability check for SHIP candidates: extract each candidate's own
genomic sequence and realign it against the same reference with bwa mem.
A candidate sitting in a repeat or segmental duplication will multi-map
(secondary/supplementary alignments, or a low MAPQ) rather than aligning
uniquely back to its own coordinates - this is a direct, self-contained
proxy for repeat content / mappability that doesn't require a separate
RepeatMasker run or external mappability track (see 05_SHIP/README.md).

Writes one row per candidate: mapq, has_secondary, has_supplementary,
n_alignments, low_mappability (mapq < --min-mapq or any secondary/
supplementary alignment).
"""
import argparse
import subprocess
import tempfile
import os
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--min-mapq", type=int, default=30)
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.candidates, encoding="utf-8") as f:
        candidates = list(csv.DictReader(f, delimiter="\t"))

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "candidates.fa")
        with open(fasta_path, "w", encoding="utf-8") as fa:
            for i, c in enumerate(candidates):
                region = f"{c['chrom']}:{int(c['start']) + 1}-{c['end']}"
                seq = subprocess.run(
                    ["samtools", "faidx", args.genome_fasta, region],
                    capture_output=True, text=True, check=True,
                ).stdout
                seq_lines = seq.splitlines()[1:]
                fa.write(f">{i}\n" + "".join(seq_lines) + "\n")

        sam_path = os.path.join(tmp, "candidates.sam")
        with open(sam_path, "w", encoding="utf-8") as sam_out:
            subprocess.run(
                ["bwa", "mem", "-t", str(args.nproc), args.genome_fasta, fasta_path],
                stdout=sam_out, stderr=subprocess.DEVNULL, check=True,
            )

        results = {}
        with open(sam_path, encoding="utf-8") as sam:
            for line in sam:
                if line.startswith("@"):
                    continue
                fields = line.split("\t")
                qname, flag, mapq = fields[0], int(fields[1]), int(fields[4])
                idx = int(qname)
                is_secondary = bool(flag & 256)
                is_supplementary = bool(flag & 2048)
                has_xa = any(f.startswith("XA:Z:") for f in fields[11:])
                r = results.setdefault(idx, {"mapq": mapq, "has_secondary": False,
                                              "has_supplementary": False, "has_xa": False,
                                              "n_alignments": 0})
                r["n_alignments"] += 1
                if not is_secondary and not is_supplementary:
                    r["mapq"] = mapq
                    r["has_xa"] = has_xa
                r["has_secondary"] = r["has_secondary"] or is_secondary
                r["has_supplementary"] = r["has_supplementary"] or is_supplementary

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("chrom\tstart\tend\tmapq\thas_secondary\thas_supplementary\thas_xa\tn_alignments\tlow_mappability\n")
        n_low = 0
        for i, c in enumerate(candidates):
            r = results.get(i, {"mapq": 0, "has_secondary": False, "has_supplementary": False,
                                 "has_xa": False, "n_alignments": 0})
            low = r["mapq"] < args.min_mapq or r["has_secondary"] or r["has_supplementary"] or r["has_xa"]
            n_low += low
            f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{r['mapq']}\t{r['has_secondary']}\t"
                     f"{r['has_supplementary']}\t{r['has_xa']}\t{r['n_alignments']}\t{low}\n")

    print(f"Checked {len(candidates)} candidates, {n_low} flagged low-mappability -> {args.out}")


if __name__ == "__main__":
    main()
