#!/usr/bin/env python3
"""
Build a canine-usable risk-gene list (oncogenes, tumor suppressors, core
essential genes) for the 05_SHIP hard-veto step, from two open, citable
sources:

  - CancerMine (Lever et al. 2019, Nature Methods; CC0) - a literature-mined
    database of gene roles in cancer across cancer types, with a citation
    count per (gene, role, cancer-type) triple. We take the union of
    citation counts per (gene, role) across all cancer types and keep only
    genes at or above --min-citations, to filter out single-mention/
    low-confidence text-mining hits.
  - CEG2 / CEGv2 (Hart et al. 2017, G3; via github.com/hart-lab/bagel) -
    684 human core essential genes from genome-wide CRISPR screens across
    multiple cell lines.

Both lists use human HGNC gene symbols. Since named genes in the
ROS_Cfam_1.0 RefSeq annotation already carry the same symbol as their human
ortholog (e.g. TSC22D3, PRPS1, RIT2 - see 05_SHIP/ship_raw_candidates.tsv),
a direct symbol match against SHIP's flanking-gene names is a reasonable
first-pass proxy for canine orthology. This is a starter list, not an
exhaustive canine-curated one - see 05_SHIP/README.md "Known gaps".
"""
import argparse
import csv
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ceg2", required=True, help="CEGv2.txt (Hart et al. 2017)")
    ap.add_argument("--cancermine", required=True, help="cancermine_collated.tsv")
    ap.add_argument("--min-citations", type=int, default=5,
                     help="minimum summed citation count across cancer types to keep a gene/role")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    genes = {}  # symbol -> set of categories

    with open(args.ceg2, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genes.setdefault(row["GENE"], set()).add("essential")

    citations = defaultdict(int)  # (gene, role) -> summed citation_count
    with open(args.cancermine, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["role"] not in ("Oncogene", "Tumor_Suppressor"):
                continue
            citations[(row["gene_normalized"], row["role"])] += int(row["citation_count"])

    for (gene, role), total in citations.items():
        if total >= args.min_citations:
            genes.setdefault(gene, set()).add(role.lower())

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("gene_symbol\tcategories\n")
        for gene in sorted(genes):
            f.write(f"{gene}\t{','.join(sorted(genes[gene]))}\n")

    n_essential = sum(1 for c in genes.values() if "essential" in c)
    n_onco = sum(1 for c in genes.values() if "oncogene" in c)
    n_tsg = sum(1 for c in genes.values() if "tumor_suppressor" in c)
    print(f"Wrote {len(genes)} risk genes to {args.out}")
    print(f"  essential: {n_essential}, oncogene: {n_onco}, tumor_suppressor: {n_tsg}")
    print(f"  (min_citations={args.min_citations} for oncogene/tumor_suppressor)")


if __name__ == "__main__":
    main()
