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
        c["hard_veto"] = (c["veto_tad_boundary"] or c["veto_atac_peak"]
                           or c["veto_risk_gene"] or c["veto_low_mappability"])
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
                  "veto_low_mappability", "no_rrbs_coverage",
                  "atac_mean", "atac_variability", "rrbs_mean", "rrbs_variability",
                  "tad_boundary_distance",
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
          f"{sum(1 for c in candidates if c['veto_low_mappability'])} low mappability), "
          f"{len(survivors)} ranked and passing.")
    print(f"Wrote full table to {args.out_scored}")
    print(f"Wrote ranked passing candidates BED to {args.out_passing_bed}")
    if survivors:
        top = survivors[0]
        print(f"Top candidate: {top['chrom']}:{top['start']}-{top['end']} "
              f"({top['orientation']}, score={top['final_score']})")


if __name__ == "__main__":
    main()
