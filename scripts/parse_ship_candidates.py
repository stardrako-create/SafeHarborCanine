#!/usr/bin/env python3
"""
Parse a SHIP (Leitao et al. 2025, github.com/MCLeitao/Ship) result .txt file
into a clean candidates TSV/BED. SHIP's own output interleaves, per
candidate: two Python-repr'd GFF rows for the flanking genes, a separator, a
"Length = Nbp" line, a FASTA header, and the intergenic sequence itself -
not something downstream tools can consume directly.

The candidate interval is the gap between the two flanking genes (left
gene's end+1 to right gene's start-1). Orientation is inferred from the
flanking strands, matching SHIP's own convergent/divergent/tandem
convention: (+,-) -> convergent, (-,+) -> divergent, same-strand -> tandem.
"""
import argparse
import ast
import re


def parse_gene_line(line):
    row = ast.literal_eval(line.strip())
    chrom, _source, _feature, start, end, _score, strand, _phase, attrs, _code = row
    gene_id_match = re.search(r"ID=gene-([^;]+)", attrs)
    gene_id = gene_id_match.group(1) if gene_id_match else "NA"
    return {"chrom": chrom, "start": int(start), "end": int(end), "strand": strand, "gene_id": gene_id}


def orientation(left_strand, right_strand):
    if left_strand == "+" and right_strand == "-":
        return "convergent"
    if left_strand == "-" and right_strand == "+":
        return "divergent"
    return "tandem"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ship-result", required=True, help="SHIP *_result_SHIP.txt file")
    ap.add_argument("--out", required=True, help="output TSV path")
    args = ap.parse_args()

    candidates = []
    pending_genes = []
    with open(args.ship_result, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("['") or line.startswith('["'):
                pending_genes.append(parse_gene_line(line))
            elif line.startswith("Length = "):
                if len(pending_genes) != 2:
                    raise SystemExit(f"Expected 2 flanking genes before a Length line, got {len(pending_genes)}")
                g1, g2 = sorted(pending_genes, key=lambda g: g["start"])
                candidates.append({
                    "chrom": g1["chrom"],
                    "start": g1["end"],  # 0-based BED start = 1-based end of left gene
                    "end": g2["start"] - 1,
                    "orientation": orientation(g1["strand"], g2["strand"]),
                    "left_gene": g1["gene_id"],
                    "left_strand": g1["strand"],
                    "right_gene": g2["gene_id"],
                    "right_strand": g2["strand"],
                })
                pending_genes = []

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("chrom\tstart\tend\tlength\torientation\tleft_gene\tleft_strand\tright_gene\tright_strand\n")
        for c in candidates:
            length = c["end"] - c["start"]
            f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{length}\t{c['orientation']}\t"
                    f"{c['left_gene']}\t{c['left_strand']}\t{c['right_gene']}\t{c['right_strand']}\n")

    print(f"Parsed {len(candidates)} SHIP candidates -> {args.out}")
    orientations = {}
    for c in candidates:
        orientations[c["orientation"]] = orientations.get(c["orientation"], 0) + 1
    print(f"Orientation breakdown: {orientations}")


if __name__ == "__main__":
    main()
