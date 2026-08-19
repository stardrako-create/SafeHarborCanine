#!/usr/bin/env python3
"""
Filter and rank SHIP-generated intergenic candidates (parse_ship_candidates.py
output) against our own tracks, instead of SHIP's own UCSC/Ensembl/regulatory
cross-reference (built for model organisms, not adapted for dog - see
scripts/build_hic_tracks.py header and the project README for why).

Two-stage design mirrors Leitao et al. 2025 (SHIP, candidate generation) +
Shrestha et al. 2022 (GEG-SH, tissue-specific biological filtering):

  hard vetoes (exclude outright)
    - overlaps a Hi-C TAD boundary (04_tracks_processadas/.../HiC/tad_boundaries.bed)
    - overlaps an ATAC consensus peak (.../ATAC/mother_track/consensus_peaks_ATAC.bed) -
      an unannotated regulatory element sitting inside a nominally intergenic window
    - either flanking gene is a known oncogene, tumor suppressor, or core
      essential gene (scripts/build_canine_risk_genes.py; CancerMine +
      Hart et al. 2017 CEG2, symbol-matched as a starter proxy for canine
      orthology - see 05_SHIP/README.md "Known gaps")
    - self-mappability check flags the candidate (scripts/check_candidate_mappability.py:
      the candidate's own sequence realigned against the genome comes back
      ambiguous - MAPQ < 30, or a secondary/supplementary/XA multi-mapping
      hit - a proxy for repeat content without a full RepeatMasker run)
    - within --mirna-min-distance (default 300kb) of an annotated miRNA
      (scripts/extract_gff3_features.py --feature-types miRNA; Ahmed et al.
      2026, Cells, criterion 3 - see 05_SHIP/ahmed2026_checklist_comparison.md)
    - within --cancer-gene-radius (default 300kb) of ANY gene in the risk-gene
      list, not just the two SHIP flanking genes (criterion 2 - a genome-wide
      radius search, using canine_all_genes.bed against canine_risk_genes.tsv)
    - a THIRD gene (beyond the two that define the SHIP window) sits within
      --gene_dense_radius (default 50kb) of either window edge (criterion 1 -
      note: SHIP windows are 50-75kb and touch genes at both edges by
      construction, so a literal "any position >=50kb from every gene" is
      structurally unsatisfiable for a window this size; this instead checks
      whether a THIRD gene crowds the neighborhood beyond the two that
      already define the window, which is the real safety concern the
      criterion is protecting against - see 05_SHIP/ahmed2026_checklist_comparison.md)
    - overlaps an annotated lncRNA or small RNA gene (lnc_RNA, tRNA, snoRNA,
      snRNA, guide_RNA, rRNA, SRP_RNA, RNase_P_RNA - criterion 6; miRNA is
      handled separately above as a 300kb radius, criterion 3)
    - the candidate's own TAD (scripts/build_tad_intervals.py - the span
      between two consecutive merged boundary calls) contains ANY risk gene
      anywhere in the domain, not just the two flanking genes or genes
      within a fixed radius (criterion 8's actual intent: a 3D-organization
      check, broader than the linear-distance checks above)

  soft score (rank what's left), each component min-max normalized to [0,1]
  across the surviving candidates, then averaged with user-settable weights:
    - stability_atac    = 1 - norm(ATAC variability across 71 dogs)
    - stability_rrbs     = 1 - norm(RRBS variability across 71 dogs)
    - low_methylation    = 1 - norm(RRBS weighted-mean methylation)
    - tad_distance       =     norm(distance to nearest TAD boundary)
    - moderate_atac      = 1 - 2*|norm(ATAC weighted mean) - 0.5|
                            (peaks at the population median accessibility -
                            SHIP already picked intergenic windows, so neither
                            fully closed nor unusually open is what we want)

Nothing here is a black box: every component is a plain, documented,
min-max-normalized track value: read the numbers in the output TSV directly.

Candidates with no RRBS coverage over the interval are flagged
(no_rrbs_coverage) rather than silently imputed - the RRBS-derived
components are left as NaN and excluded from that candidate's score average.
"""
import argparse
import csv

import pyBigWig


def load_bed_intervals(path):
    intervals = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields[1].isdigit():
                continue  # header row
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            intervals.setdefault(chrom, []).append((start, end))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def load_bed_with_names(path):
    intervals = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields[1].isdigit():
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name = fields[3] if len(fields) > 3 else "NA"
            intervals.setdefault(chrom, []).append((start, end, name))
    return intervals


def distance_to_nearest_named(intervals_by_chrom, chrom, start, end, name_set):
    mid = (start + end) // 2
    best = None
    for s, e, name in intervals_by_chrom.get(chrom, []):
        if name not in name_set:
            continue
        if e < mid:
            d = mid - e
        elif s > mid:
            d = s - mid
        else:
            d = 0
        if best is None or d < best:
            best = d
    return best if best is not None else float("inf")


def distance_from_point_excluding(intervals_by_chrom, chrom, point, exclude_names):
    best = None
    for s, e, name in intervals_by_chrom.get(chrom, []):
        if name in exclude_names:
            continue
        if e < point:
            d = point - e
        elif s > point:
            d = s - point
        else:
            d = 0
        if best is None or d < best:
            best = d
    return best if best is not None else float("inf")


def find_containing_interval(intervals_by_chrom, chrom, point):
    for s, e in intervals_by_chrom.get(chrom, []):
        if s <= point <= e:
            return (s, e)
    return None


def any_named_overlaps(intervals_by_chrom, chrom, start, end, name_set):
    for s, e, name in intervals_by_chrom.get(chrom, []):
        if name in name_set and s < end and e > start:
            return True
    return False


def overlaps(intervals_by_chrom, chrom, start, end):
    for s, e in intervals_by_chrom.get(chrom, []):
        if s < end and e > start:
            return True
    return False


def distance_to_nearest(intervals_by_chrom, chrom, start, end):
    mid = (start + end) // 2
    best = None
    for s, e in intervals_by_chrom.get(chrom, []):
        if e < mid:
            d = mid - e
        elif s > mid:
            d = s - mid
        else:
            d = 0
        if best is None or d < best:
            best = d
    return best if best is not None else float("inf")


def bw_mean(bw, chrom, start, end):
    try:
        val = bw.stats(chrom, start, end, type="mean")[0]
    except RuntimeError:
        return None
    return val


def minmax_norm(values):
    finite = [v for v in values if v is not None]
    if not finite:
        return [None for _ in values]
    lo, hi = min(finite), max(finite)
    if hi - lo < 1e-12:
        return [1.0 if v is not None else None for v in values]
    return [(v - lo) / (hi - lo) if v is not None else None for v in values]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="parse_ship_candidates.py output TSV")
    ap.add_argument("--atac-mean-bw", required=True)
    ap.add_argument("--atac-variability-bw", required=True)
    ap.add_argument("--atac-peaks-bed", required=True)
    ap.add_argument("--rrbs-mean-bw", required=True)
    ap.add_argument("--rrbs-variability-bw", required=True)
    ap.add_argument("--rrbs-coverage-bw", required=True)
    ap.add_argument("--tad-boundaries-bed", required=True)
    ap.add_argument("--risk-genes", required=True, help="build_canine_risk_genes.py output TSV")
    ap.add_argument("--mappability-tsv", required=True, help="check_candidate_mappability.py output TSV")
    ap.add_argument("--mirna-bed", required=True, help="extract_gff3_features.py --feature-types miRNA output")
    ap.add_argument("--mirna-min-distance", type=int, default=300_000)
    ap.add_argument("--all-genes-bed", required=True, help="extract_gff3_features.py --feature-types gene pseudogene output")
    ap.add_argument("--cancer-gene-radius", type=int, default=300_000)
    ap.add_argument("--gene-dense-radius", type=int, default=50_000)
    ap.add_argument("--lncrna-smallrna-bed", required=True,
                     help="extract_gff3_features.py lncRNA/small RNA output")
    ap.add_argument("--tad-intervals-bed", required=True, help="build_tad_intervals.py output")
    ap.add_argument("--w-stability-atac", type=float, default=1.0)
    ap.add_argument("--w-stability-rrbs", type=float, default=1.0)
    ap.add_argument("--w-low-methylation", type=float, default=1.0)
    ap.add_argument("--w-tad-distance", type=float, default=1.0)
    ap.add_argument("--w-moderate-atac", type=float, default=1.0)
    ap.add_argument("--out-scored", required=True)
    ap.add_argument("--out-passing-bed", required=True)
    args = ap.parse_args()

    with open(args.candidates, encoding="utf-8") as f:
        candidates = list(csv.DictReader(f, delimiter="\t"))
    for c in candidates:
        c["start"] = int(c["start"])
        c["end"] = int(c["end"])

    atac_peaks = load_bed_intervals(args.atac_peaks_bed)
    tad_boundaries = load_bed_intervals(args.tad_boundaries_bed)
    mirnas = load_bed_intervals(args.mirna_bed)
    all_genes = load_bed_with_names(args.all_genes_bed)
    lncrna_smallrna = load_bed_intervals(args.lncrna_smallrna_bed)
    tad_intervals = load_bed_intervals(args.tad_intervals_bed)

    risk_genes = set()
    with open(args.risk_genes, encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            risk_genes.add(line.split("\t")[0])

    low_mappability = {}
    with open(args.mappability_tsv, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            low_mappability[(row["chrom"], int(row["start"]), int(row["end"]))] = row["low_mappability"] == "True"

    atac_mean_bw = pyBigWig.open(args.atac_mean_bw)
    atac_var_bw = pyBigWig.open(args.atac_variability_bw)
    rrbs_mean_bw = pyBigWig.open(args.rrbs_mean_bw)
    rrbs_var_bw = pyBigWig.open(args.rrbs_variability_bw)
    rrbs_cov_bw = pyBigWig.open(args.rrbs_coverage_bw)

    for c in candidates:
        chrom, start, end = c["chrom"], c["start"], c["end"]
        c["veto_tad_boundary"] = overlaps(tad_boundaries, chrom, start, end)
        c["veto_atac_peak"] = overlaps(atac_peaks, chrom, start, end)
        c["veto_risk_gene"] = c["left_gene"] in risk_genes or c["right_gene"] in risk_genes
        c["veto_low_mappability"] = low_mappability.get((chrom, start, end), False)
        c["mirna_distance"] = distance_to_nearest(mirnas, chrom, start, end)
        c["veto_mirna_nearby"] = c["mirna_distance"] < args.mirna_min_distance
        c["risk_gene_radius_distance"] = distance_to_nearest_named(all_genes, chrom, start, end, risk_genes)
        c["veto_risk_gene_radius"] = c["risk_gene_radius_distance"] < args.cancer_gene_radius
        exclude = {c["left_gene"], c["right_gene"]}
        left_clear = distance_from_point_excluding(all_genes, chrom, start, exclude)
        right_clear = distance_from_point_excluding(all_genes, chrom, end, exclude)
        c["gene_dense_clearance"] = min(left_clear, right_clear)
        c["veto_gene_dense_neighborhood"] = c["gene_dense_clearance"] < args.gene_dense_radius
        c["veto_lncrna_smallrna"] = overlaps(lncrna_smallrna, chrom, start, end)
        mid = (start + end) // 2
        own_tad = find_containing_interval(tad_intervals, chrom, mid)
        if own_tad is not None:
            c["veto_tad_risk_gene"] = any_named_overlaps(all_genes, chrom, own_tad[0], own_tad[1], risk_genes)
        else:
            c["veto_tad_risk_gene"] = False  # no TAD interval could be resolved (e.g. chromosome end)
        c["hard_veto"] = (c["veto_tad_boundary"] or c["veto_atac_peak"]
                           or c["veto_risk_gene"] or c["veto_low_mappability"]
                           or c["veto_mirna_nearby"] or c["veto_risk_gene_radius"]
                           or c["veto_gene_dense_neighborhood"] or c["veto_lncrna_smallrna"]
                           or c["veto_tad_risk_gene"])
        c["tad_boundary_distance"] = distance_to_nearest(tad_boundaries, chrom, start, end)
        c["atac_mean"] = bw_mean(atac_mean_bw, chrom, start, end)
        c["atac_variability"] = bw_mean(atac_var_bw, chrom, start, end)
        c["rrbs_mean"] = bw_mean(rrbs_mean_bw, chrom, start, end)
        c["rrbs_variability"] = bw_mean(rrbs_var_bw, chrom, start, end)
        c["rrbs_coverage"] = bw_mean(rrbs_cov_bw, chrom, start, end)
        c["no_rrbs_coverage"] = not c["rrbs_coverage"]

    survivors = [c for c in candidates if not c["hard_veto"]]

    atac_var_norm = minmax_norm([c["atac_variability"] for c in survivors])
    atac_mean_norm = minmax_norm([c["atac_mean"] for c in survivors])
    rrbs_mean_norm = minmax_norm([c["rrbs_mean"] if not c["no_rrbs_coverage"] else None for c in survivors])
    rrbs_var_norm = minmax_norm([c["rrbs_variability"] if not c["no_rrbs_coverage"] else None for c in survivors])
    tad_dist_norm = minmax_norm([c["tad_boundary_distance"] for c in survivors])

    weights = {
        "stability_atac": args.w_stability_atac,
        "stability_rrbs": args.w_stability_rrbs,
        "low_methylation": args.w_low_methylation,
        "tad_distance": args.w_tad_distance,
        "moderate_atac": args.w_moderate_atac,
    }

    for c, av, mv, rm, rv, td in zip(survivors, atac_var_norm, atac_mean_norm, rrbs_mean_norm, rrbs_var_norm, tad_dist_norm):
        components = {}
        components["stability_atac"] = (1.0 - av) if av is not None else None
        components["moderate_atac"] = (1.0 - 2.0 * abs(mv - 0.5)) if mv is not None else None
        components["low_methylation"] = (1.0 - rm) if rm is not None else None
        components["stability_rrbs"] = (1.0 - rv) if rv is not None else None
        components["tad_distance"] = td if td is not None else None

        weighted_sum, weight_total = 0.0, 0.0
        for name, value in components.items():
            if value is None:
                continue
            weighted_sum += weights[name] * value
            weight_total += weights[name]
            c[f"score_{name}"] = round(value, 4)
        c["final_score"] = round(weighted_sum / weight_total, 4) if weight_total > 0 else None

    survivors.sort(key=lambda c: (c["final_score"] is None, -(c["final_score"] or 0)))

    fieldnames = ["chrom", "start", "end", "length", "orientation", "left_gene", "right_gene",
                  "hard_veto", "veto_tad_boundary", "veto_atac_peak", "veto_risk_gene",
                  "veto_low_mappability", "veto_mirna_nearby", "veto_risk_gene_radius",
                  "veto_gene_dense_neighborhood", "veto_lncrna_smallrna", "veto_tad_risk_gene", "no_rrbs_coverage",
                  "atac_mean", "atac_variability", "rrbs_mean", "rrbs_variability",
                  "tad_boundary_distance", "mirna_distance", "risk_gene_radius_distance", "gene_dense_clearance",
                  "score_stability_atac", "score_moderate_atac", "score_low_methylation",
                  "score_stability_rrbs", "score_tad_distance", "final_score"]

    with open(args.out_scored, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for c in candidates:
            writer.writerow(c)

    with open(args.out_passing_bed, "w", encoding="utf-8") as f:
        for c in survivors:
            name = f"{c['orientation']}_{c['left_gene']}_{c['right_gene']}_score{c['final_score']}"
            f.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{name}\t{c['final_score']}\n")

    n_veto = sum(1 for c in candidates if c["hard_veto"])
    print(f"{len(candidates)} candidates total, {n_veto} hard-vetoed "
          f"({sum(1 for c in candidates if c['veto_tad_boundary'])} TAD boundary, "
          f"{sum(1 for c in candidates if c['veto_atac_peak'])} ATAC peak overlap, "
          f"{sum(1 for c in candidates if c['veto_risk_gene'])} risk gene, "
          f"{sum(1 for c in candidates if c['veto_low_mappability'])} low mappability, "
          f"{sum(1 for c in candidates if c['veto_mirna_nearby'])} miRNA nearby, "
          f"{sum(1 for c in candidates if c['veto_risk_gene_radius'])} risk gene within radius, "
          f"{sum(1 for c in candidates if c['veto_gene_dense_neighborhood'])} gene-dense neighborhood, "
          f"{sum(1 for c in candidates if c['veto_lncrna_smallrna'])} lncRNA/smallRNA overlap, "
          f"{sum(1 for c in candidates if c['veto_tad_risk_gene'])} risk gene in own TAD), "
          f"{len(survivors)} ranked and passing.")
    print(f"Wrote full table to {args.out_scored}")
    print(f"Wrote ranked passing candidates BED to {args.out_passing_bed}")
    if survivors:
        top = survivors[0]
        print(f"Top candidate: {top['chrom']}:{top['start']}-{top['end']} "
              f"({top['orientation']}, score={top['final_score']})")


if __name__ == "__main__":
    main()
