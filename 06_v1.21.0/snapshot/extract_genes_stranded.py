#!/usr/bin/env python3
"""
Extract protein-coding gene intervals WITH strand from the GFF3, written
2026-09-11 to fix a real literature-fidelity bug: Ahmed et al. 2026's
criterion 1 is explicit in its body text (Section 2, page 3) - "SHSs must
be located at a safe distance, typically at least 50 kilobases (kb), from
the 5' end of coding genes" - not from the gene body generally. The
existing canine_all_genes.bed has no strand column, so score_ship_
candidates_v2.py's gene-density check was measuring distance to the
nearest gene-body edge instead, which is a different criterion, not just a
conservative version of the same one (confirmed by re-reading the paper's
actual PDF text directly, not the project's own paraphrased summary).

Filtered to gene_biotype=protein_coding specifically, matching Ahmed's
"coding genes" wording - miRNA and lncRNA/smallRNA genes are already
covered by their own separate criteria (3 and 6) elsewhere in this
pipeline; including them here would double up on those, not implement
criterion 1.

Output: chrom, start(0-based), end, name, score(0, unused), strand - a
proper 6-column BED.
"""
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args()

    n_total_genes = 0
    n_written = 0
    with open(args.gff3, encoding="utf-8") as fin, open(args.out_bed, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            n_total_genes += 1
            attrs = fields[8]
            if "gene_biotype=protein_coding" not in attrs:
                continue
            chrom = fields[0]
            start_1based, end = int(fields[3]), int(fields[4])
            strand = fields[6]
            if strand not in ("+", "-"):
                continue  # skip anything without a definite strand - can't compute a 5' end
            name = None
            for part in attrs.split(";"):
                if part.startswith("Name="):
                    name = part[len("Name="):]
                    break
            if name is None:
                continue
            start_0based = start_1based - 1
            fout.write(f"{chrom}\t{start_0based}\t{end}\t{name}\t0\t{strand}\n")
            n_written += 1

    print(f"{n_total_genes} total 'gene' features in GFF3, "
          f"{n_written} written as protein-coding, stranded, named genes to {args.out_bed}")


if __name__ == "__main__":
    main()
