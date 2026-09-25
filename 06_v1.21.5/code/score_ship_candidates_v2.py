#!/usr/bin/env python3
"""
v1.20.0 scoring - identical logic to score_ship_candidates.py (V12, current),
pointed at the unified 06_v1.20.0/ track locations instead of the old
04_tracks_processadas/.../ATAC*/RRBS split-cohort layout, and using
bw_utils.py's shared implicit-zero handling instead of a locally
reimplemented copy (see bw_utils.py's own docstring for why that matters -
that duplication is exactly what caused the 2026-09-08 peak_frequency
background bug in the first place).

Every hard veto, the soft-score formula, and the final_score_percentile
display framing are unchanged from V12 - this file exists to point at
v1.20.0's tracks and consolidate the shared utility, not to change the
methodology. See score_ship_candidates.py's own docstring for the full,
still-accurate description of every veto and score component.

**One real fix, 2026-09-10**: the ATAC accessibility-floor veto
(veto_low_atac_accessibility) was comparing a window-averaged candidate
value (mean over the full 50-75kb SHIP span) against a background built
from individual 25bp bins - a statistical mismatch (window-averaging
mechanically compresses variance vs. individual bins), which made nearly
every candidate look accessible-enough regardless of true position:
254/461 failed under the old (buggy) comparison. Isolated with a
4-condition test - narrowing the candidate's own window alone didn't fix
it (269/461 still failed); matching the background's statistic to the
candidate's did (9/461 failed) - the bug was the mismatch, not the window
width. Fixed properly: the candidate must now pass TWO separate,
internally-consistent checks - build_matched_window_background()'s wide
check (neighborhood-level openness) and build_narrow_window_background()'s
narrow check (+/-10kb around the window center, a proxy for the likely
eventual insertion point, not yet chosen pending Vasco's vector decision).
Each compares like-for-like; neither is cross-compared with the other's
background. Result: 13/461 fail, 26 candidates pass every criterion
(recovering the same set V11 had before the original, buggy floor was
added - including NC_051807.1:10779807-10837318, the cross-validated
convergence point with Ehsan's original 27).

**Same bug, remaining 4 components, 2026-09-10 (cycle 2)**: a code review
found the accessibility fix above had not been extended to
atac_variability, rrbs_mean, rrbs_variability, or atac_peak_frequency -
each was still a window statistic compared against a background of
individual bins (bw_utils.build_bw_background /
build_sparse_bg_with_implicit_zero, now removed from this file).
atac_peak_frequency had a second, compounding problem: it's a MAX over
~2000-3000 bins per candidate window compared against single-bin draws - a
window MAX is a biased-high order statistic versus individual draws by
construction, independent of the mismatch bug. Fixed the same way as
accessibility: every candidate-side value and its background now go
through bw_utils.track_value()/build_matched_window_background(), so the
statistic can no longer silently diverge between the two sides. tad_distance
was checked and is NOT affected - both its candidate value and its
background were already the same statistic (distance from a point to the
nearest TAD boundary). Also fixed in the same pass: build_matched_window_
background()'s chromosome sampling was uniform-per-chromosome rather than
length-weighted, over-representing small contigs; and every bw.stats() call
in this file now passes exact=True, so pyBigWig never silently substitutes
a coarser zoom-level summary for the exact interval statistic.

**Literature-fidelity fix, 2026-09-11**: veto_gene_dense_neighborhood did
NOT implement Ahmed et al. 2026's actual criterion 1 - verified directly
against the paper's PDF (Section 2, page 3): "SHSs must be located at a
safe distance, typically at least 50 kilobases (kb), from the 5' end of
coding genes." The prior implementation measured distance to the nearest
gene-BODY edge, excluding the two flanking genes by name - a different
criterion, not a conservative version of the same one. The exclusion was
justified in 05_SHIP/ahmed2026_checklist_comparison.md by an argument that
a 50-75kb window touching genes at both edges can never be 50kb from
either flanking gene - true for body-distance, but SHIP only returns
convergent gene pairs, so the flanking genes' near edges (facing the
window) are their 3' ends by construction; their 5' ends face outward and
aren't automatically far enough just because they happen to define the
window. Fixed: distance_to_nearest_5prime_end() measures to the nearest
5' end (strand-aware, via extract_genes_stranded.py's protein-coding-only
gene set) among ALL genes, no exclusion. The superseded body-distance
metric is kept as gene_dense_clearance_old_bodydist for comparison, not
used in hard_veto. Separately, veto_external_regulatory_element and
veto_atac_peak/low_peak_frequency are project-specific additions beyond
Ahmed's 8 literal criteria, not implementations of any of them - say
"passes every implemented veto," not "passes every Ahmed criterion."

**RRBS evidence-masking fix, 2026-09-12 (found by Codex, a separate working
session on this same project - see 06_v1.21.1/, adopted here)**: rrbs_mean
was a plain window average over methylation_weighted_mean.bw, which zero-
pads bins with no dog coverage (dense-track convention) - the same
implicit-zero-vs-missing conflation as bug 1 above, just at the read level
instead of the track-build level. A 50-75kb window can be ~88% uncovered
bins (RRBS is sparse, concentrated at MspI sites), so this silently
diluted rrbs_mean toward 0 regardless of the true methylation level at
positions that actually had data. Fixed via bw_utils.evidence_summary():
the mean is now taken only over bases the coverage track confirms were
observed, on both the candidate side and (build_matched_window_
background()'s coverage_bw param) the background side. rrbs_observed_
fraction/rrbs_observed_bp are new columns exposing how much of each
candidate's window this mean actually rests on - low values mean "sparse
evidence," not necessarily "unusual methylation."
"""
import argparse
import csv
import math

import numpy as np
import pyBigWig

from bw_utils import (track_value, build_matched_window_background,
                       build_narrow_window_background, evidence_summary)


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


def load_genes_stranded(path):
    """chrom -> [(start, end, name, strand, five_prime_end), ...] - from
    extract_genes_stranded.py's output (protein-coding only, real strand).
    five_prime_end is start for '+' genes, end for '-' genes - this is what
    Ahmed et al. 2026's criterion 1 actually measures distance to (verified
    against the paper directly, 2026-09-11), not the nearest gene-body edge."""
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
            strand = fields[5] if len(fields) > 5 else "+"
            five_prime = start if strand == "+" else end
            intervals.setdefault(chrom, []).append((start, end, name, strand, five_prime))
    return intervals


def distance_to_nearest_5prime_end(genes_by_chrom, chrom, point):
    """Distance from a point to the nearest coding gene's 5' end, among
    ALL genes - not excluding the two flanking genes. SHIP only returns
    convergent gene pairs, so the flanking genes' NEAR edges (facing the
    candidate window) are their 3' ends by construction; their 5' ends
    face outward and are not guaranteed to already clear 50kb just because
    they define the window - if a flanking gene is short enough that its
    own 5' end is still close, Ahmed's literal criterion doesn't exempt it
    just because it happens to be one of the two genes SHIP anchored on."""
    best = None
    for s, e, name, strand, five_prime in genes_by_chrom.get(chrom, []):
        d = abs(point - five_prime)
        if best is None or d < best:
            best = d
    return best if best is not None else float("inf")


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


def narrow_window_mean(bw, chrom, start, end, radius=10_000):
    """The candidate's own accessibility right around its likely eventual
    insertion point (window center) - a 20kb window (+/-10kb), not the
    full 50-75kb SHIP span. The exact insertion coordinate isn't chosen
    yet (depends on Vasco's vector/backbone choice), so this is a
    reasonable proxy centered on the window, not a final answer."""
    mid = (start + end) // 2
    return track_value(bw, chrom, max(0, mid - radius), mid + radius, stat_type="mean")


def build_tad_distance_background(tad_boundaries, chrom_sizes, step=1000):
    chunks = []
    for chrom, size in chrom_sizes.items():
        grid = np.arange(0, size, step, dtype=np.int64)
        boundaries = tad_boundaries.get(chrom, [])
        if not boundaries:
            chunks.append(np.full(grid.shape, size, dtype=np.float64))
            continue
        starts = np.array([s for s, _ in boundaries], dtype=np.int64)
        ends = np.array([e for _, e in boundaries], dtype=np.int64)
        idx = np.searchsorted(starts, grid, side="right") - 1
        idx = np.clip(idx, 0, len(starts) - 1)
        idx_next = np.clip(idx + 1, 0, len(starts) - 1)
        d_left = np.maximum(0, np.maximum(starts[idx] - grid, grid - ends[idx]))
        d_right = np.maximum(0, np.maximum(starts[idx_next] - grid, grid - ends[idx_next]))
        chunks.append(np.minimum(d_left, d_right).astype(np.float64))
    background = np.concatenate(chunks) if chunks else np.array([])
    background.sort()
    return background


def percentile_rank(sorted_background, value):
    if value is None or sorted_background.size == 0:
        return None
    if isinstance(value, float) and math.isinf(value):
        return 1.0
    idx = np.searchsorted(sorted_background, value, side="right")
    return float(idx) / sorted_background.size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--atac-mean-bw", required=True)
    ap.add_argument("--atac-variability-bw", required=True)
    ap.add_argument("--atac-peaks-bed", required=True)
    ap.add_argument("--atac-peak-frequency-bw", required=True)
    ap.add_argument("--min-atac-accessibility-percentile", type=float, default=0.55)
    ap.add_argument("--rrbs-mean-bw", required=True)
    ap.add_argument("--rrbs-variability-bw", required=True)
    ap.add_argument("--rrbs-coverage-bw", required=True)
    ap.add_argument("--tad-boundaries-bed", required=True)
    ap.add_argument("--risk-genes", required=True)
    ap.add_argument("--mappability-tsv", required=True)
    ap.add_argument("--mirna-bed", required=True)
    ap.add_argument("--mirna-min-distance", type=int, default=300_000)
    ap.add_argument("--all-genes-bed", required=True)
    ap.add_argument("--all-genes-stranded-bed", required=True,
                     help="protein-coding genes with strand (extract_genes_stranded.py) "
                          "- used for the 5'-end-based gene-density check, since Ahmed "
                          "et al. 2026 criterion 1 is literally about distance from the "
                          "5' end of coding genes, not gene-body distance")
    ap.add_argument("--cancer-gene-radius", type=int, default=300_000)
    ap.add_argument("--gene-dense-radius", type=int, default=50_000,
                     help="Ahmed et al. 2026 criterion 1's literal value - do not "
                          "change without flagging the deviation explicitly")
    ap.add_argument("--lncrna-smallrna-bed", required=True)
    ap.add_argument("--tad-intervals-bed", required=True)
    ap.add_argument("--repeat-content-tsv", default=None)
    ap.add_argument("--repeat-content-threshold", type=float, default=50.0)
    ap.add_argument("--ultraconserved-tsv", default=None)
    ap.add_argument("--ultraconserved-threshold", type=float, default=6.5)
    ap.add_argument("--external-regulatory-bed", default=None)
    ap.add_argument("--w-stability-atac", type=float, default=1.0)
    ap.add_argument("--w-stability-rrbs", type=float, default=1.0)
    ap.add_argument("--w-low-methylation", type=float, default=1.0)
    ap.add_argument("--w-tad-distance", type=float, default=1.0)
    ap.add_argument("--w-low-peak-frequency", type=float, default=1.0)
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
    genes_stranded = load_genes_stranded(args.all_genes_stranded_bed)
    lncrna_smallrna = load_bed_intervals(args.lncrna_smallrna_bed)
    tad_intervals = load_bed_intervals(args.tad_intervals_bed)

    risk_genes = set()
    with open(args.risk_genes, encoding="utf-8") as f:
        next(f)
        for line in f:
            risk_genes.add(line.split("\t")[0])

    low_mappability = {}
    with open(args.mappability_tsv, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            low_mappability[(row["chrom"], int(row["start"]), int(row["end"]))] = row["low_mappability"] == "True"

    repeat_content = {}
    if args.repeat_content_tsv:
        with open(args.repeat_content_tsv, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                repeat_content[(row["chrom"], int(row["start"]), int(row["end"]))] = float(row["pct_repeat"])

    ultraconserved = {}
    if args.ultraconserved_tsv:
        with open(args.ultraconserved_tsv, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                ultraconserved[(row["chrom"], int(row["start"]), int(row["end"]))] = float(row["max_50bp_rolling_phyloP"])

    external_regulatory = load_bed_intervals(args.external_regulatory_bed) if args.external_regulatory_bed else {}

    atac_mean_bw = pyBigWig.open(args.atac_mean_bw)
    atac_var_bw = pyBigWig.open(args.atac_variability_bw)
    atac_freq_bw = pyBigWig.open(args.atac_peak_frequency_bw)
    rrbs_mean_bw = pyBigWig.open(args.rrbs_mean_bw)
    rrbs_var_bw = pyBigWig.open(args.rrbs_variability_bw)
    rrbs_cov_bw = pyBigWig.open(args.rrbs_coverage_bw)

    print("Building ATAC accessibility backgrounds (dual check: narrow-window + "
          "wide-window, each matched to its own candidate-side statistic)...")
    candidate_lengths = [c["end"] - c["start"] for c in candidates]
    atac_mean_bg_wide = build_matched_window_background(atac_mean_bw, candidate_lengths)
    atac_mean_bg_narrow = build_narrow_window_background(atac_mean_bw)
    print(f"  wide background: n={atac_mean_bg_wide.size}, p55={np.percentile(atac_mean_bg_wide,55):.4g}")
    print(f"  narrow background: n={atac_mean_bg_narrow.size}, p55={np.percentile(atac_mean_bg_narrow,55):.4g}")

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
        # Ahmed et al. 2026 criterion 1, literal text (verified directly
        # against the paper, Section 2, 2026-09-11): ">=50kb from the 5'
        # end of coding genes" - distance to the nearest 5' end, among ALL
        # protein-coding genes, not the nearest gene-body edge, and not
        # excluding the two flanking genes (see distance_to_nearest_5prime_
        # end()'s docstring for why the old exclusion was itself based on a
        # mistaken premise). Checked from both window edges, worst case
        # (min) - so every point inside the window clears the radius from
        # every gene's 5' end, not just the one eventually chosen as the
        # insertion site.
        left_5p = distance_to_nearest_5prime_end(genes_stranded, chrom, start)
        right_5p = distance_to_nearest_5prime_end(genes_stranded, chrom, end)
        c["gene_5prime_clearance"] = min(left_5p, right_5p)
        c["veto_gene_dense_neighborhood"] = c["gene_5prime_clearance"] < args.gene_dense_radius
        # Superseded gene-BODY-distance metric (the pre-2026-09-11 veto
        # logic) - kept only so diff_scored_candidates.py can show exactly
        # which candidates flip and why, not used in hard_veto anymore.
        exclude = {c["left_gene"], c["right_gene"]}
        left_clear = distance_from_point_excluding(all_genes, chrom, start, exclude)
        right_clear = distance_from_point_excluding(all_genes, chrom, end, exclude)
        c["gene_dense_clearance_old_bodydist"] = min(left_clear, right_clear)
        c["veto_lncrna_smallrna"] = overlaps(lncrna_smallrna, chrom, start, end)
        c["pct_repeat"] = repeat_content.get((chrom, start, end))
        c["veto_high_repeat_content"] = (
            c["pct_repeat"] is not None and c["pct_repeat"] > args.repeat_content_threshold
        )
        c["max_50bp_rolling_phyloP"] = ultraconserved.get((chrom, start, end))
        c["veto_ultraconserved_element"] = (
            c["max_50bp_rolling_phyloP"] is not None
            and c["max_50bp_rolling_phyloP"] > args.ultraconserved_threshold
        )
        c["veto_external_regulatory_element"] = (
            bool(args.external_regulatory_bed) and overlaps(external_regulatory, chrom, start, end)
        )
        c["atac_mean"] = track_value(atac_mean_bw, chrom, start, end)
        c["atac_mean_percentile"] = percentile_rank(atac_mean_bg_wide, c["atac_mean"])
        c["atac_mean_narrow"] = narrow_window_mean(atac_mean_bw, chrom, start, end)
        c["atac_mean_narrow_percentile"] = percentile_rank(atac_mean_bg_narrow, c["atac_mean_narrow"])
        fails_wide = (
            c["atac_mean_percentile"] is None
            or c["atac_mean_percentile"] < args.min_atac_accessibility_percentile
        )
        fails_narrow = (
            c["atac_mean_narrow_percentile"] is None
            or c["atac_mean_narrow_percentile"] < args.min_atac_accessibility_percentile
        )
        # Must pass BOTH: wide (overall neighborhood openness) and narrow
        # (likely-insertion-point openness) - each checked against its own
        # matched-statistic background, never cross-compared (see
        # build_matched_window_background's docstring for the bug this fixes).
        c["veto_low_atac_accessibility"] = fails_wide or fails_narrow
        mid = (start + end) // 2
        own_tad = find_containing_interval(tad_intervals, chrom, mid)
        if own_tad is not None:
            c["veto_tad_risk_gene"] = any_named_overlaps(all_genes, chrom, own_tad[0], own_tad[1], risk_genes)
        else:
            c["veto_tad_risk_gene"] = False
        c["hard_veto"] = (c["veto_tad_boundary"] or c["veto_atac_peak"]
                           or c["veto_risk_gene"] or c["veto_low_mappability"]
                           or c["veto_mirna_nearby"] or c["veto_risk_gene_radius"]
                           or c["veto_gene_dense_neighborhood"] or c["veto_lncrna_smallrna"]
                           or c["veto_tad_risk_gene"] or c["veto_high_repeat_content"]
                           or c["veto_ultraconserved_element"] or c["veto_external_regulatory_element"]
                           or c["veto_low_atac_accessibility"])
        c["tad_boundary_distance"] = distance_to_nearest(tad_boundaries, chrom, start, end)
        c["atac_variability"] = track_value(atac_var_bw, chrom, start, end)
        c["atac_peak_frequency"] = track_value(atac_freq_bw, chrom, start, end,
                                                stat_type="max", implicit_zero=True)
        c["rrbs_mean"], c["rrbs_observed_fraction"], c["rrbs_observed_bp"] = (
            evidence_summary(rrbs_mean_bw, rrbs_cov_bw, chrom, start, end)
        )
        c["rrbs_variability"] = track_value(rrbs_var_bw, chrom, start, end)
        c["rrbs_coverage"] = track_value(rrbs_cov_bw, chrom, start, end)
        c["no_rrbs_coverage"] = c["rrbs_mean"] is None

    # survivors (the "passes_recorded_checks" definition, stricter than just
    # "not hard_veto") is computed later, once missing-evidence is known.

    print("Building matched-window background distributions for percentile scoring "
          "(same statistic and window-length distribution as each component's "
          "candidate-side value - see track_value()/build_matched_window_background() "
          "in bw_utils.py for why this has to match)...")
    atac_var_bg = build_matched_window_background(atac_var_bw, candidate_lengths)
    atac_freq_bg = build_matched_window_background(atac_freq_bw, candidate_lengths,
                                                     stat_type="max", implicit_zero=True)
    rrbs_mean_bg = build_matched_window_background(rrbs_mean_bw, candidate_lengths,
                                                     coverage_bw=rrbs_cov_bw)
    rrbs_var_bg = build_matched_window_background(rrbs_var_bw, candidate_lengths)
    tad_dist_bg = build_tad_distance_background(tad_boundaries, atac_mean_bw.chroms())
    print(f"  atac_variability background: n={atac_var_bg.size}, p50={np.percentile(atac_var_bg,50):.4g}")
    print(f"  atac_peak_frequency background: n={atac_freq_bg.size}, p50={np.percentile(atac_freq_bg,50):.4g}")
    print(f"  rrbs_mean background: n={rrbs_mean_bg.size}, p50={np.percentile(rrbs_mean_bg,50):.4g}")
    print(f"  rrbs_variability background: n={rrbs_var_bg.size}, p50={np.percentile(rrbs_var_bg,50):.4g}")
    print(f"  TAD-distance background: {tad_dist_bg.size:,} grid points (unaffected by this fix - "
          f"already a matched point-to-boundary-distance statistic on both sides)")

    weights = {
        "stability_atac": args.w_stability_atac,
        "stability_rrbs": args.w_stability_rrbs,
        "low_methylation": args.w_low_methylation,
        "tad_distance": args.w_tad_distance,
        "low_peak_frequency": args.w_low_peak_frequency,
    }

    for c in candidates:
        av = percentile_rank(atac_var_bg, c["atac_variability"])
        pf = percentile_rank(atac_freq_bg, c["atac_peak_frequency"])
        rm = percentile_rank(rrbs_mean_bg, c["rrbs_mean"]) if not c["no_rrbs_coverage"] else None
        rv = percentile_rank(rrbs_var_bg, c["rrbs_variability"]) if not c["no_rrbs_coverage"] else None
        td = percentile_rank(tad_dist_bg, c["tad_boundary_distance"])

        components = {}
        components["stability_atac"] = (1.0 - av) if av is not None else None
        components["low_peak_frequency"] = (1.0 - pf) if pf is not None else None
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

    all_final_scores = np.array(
        [c["final_score"] for c in candidates if c["final_score"] is not None]
    )
    all_final_scores.sort()
    for c in candidates:
        if c["final_score"] is None:
            c["final_score_percentile"] = None
        else:
            c["final_score_percentile"] = round(
                percentile_rank(all_final_scores, c["final_score"]) * 100, 1
            )

    # Missing inputs are not biological failures and never count as a
    # completed check - adopted 2026-09-12 from Codex (06_v1.21.1/), a
    # separate working session on this project: a candidate with a hard
    # veto is "excluded"; one that clears every veto but is missing a
    # required annotation (mappability/repeat/conservation/TAD/external-
    # regulatory/score component) is "insufficient_evidence", not silently
    # counted as passing; only a candidate with every check both cleared
    # AND actually evaluated is "passes_recorded_checks".
    for c in candidates:
        missing = []
        locus = (c["chrom"], c["start"], c["end"])
        if locus not in low_mappability:
            missing.append("mappability")
        if c["pct_repeat"] is None:
            missing.append("repeat_content")
        if c["max_50bp_rolling_phyloP"] is None:
            missing.append("conservation")
        if find_containing_interval(tad_intervals, c["chrom"], (c["start"] + c["end"]) // 2) is None:
            missing.append("TAD_assignment")
        if not args.external_regulatory_bed:
            missing.append("external_regulatory_annotation")
        if c["final_score"] is None:
            missing.append("score_component")
        c["missing_evidence"] = ";".join(missing)
        c["evaluation_status"] = ("excluded" if c["hard_veto"] else
                                   "insufficient_evidence" if missing else "passes_recorded_checks")

    survivors = [c for c in candidates if c["evaluation_status"] == "passes_recorded_checks"]
    survivors.sort(key=lambda c: (c["final_score"] is None, -(c["final_score"] or 0)))

    fieldnames = ["evaluation_status", "missing_evidence",
                  "chrom", "start", "end", "length", "orientation", "left_gene", "right_gene",
                  "hard_veto", "veto_tad_boundary", "veto_atac_peak", "veto_risk_gene",
                  "veto_low_mappability", "veto_mirna_nearby", "veto_risk_gene_radius",
                  "veto_gene_dense_neighborhood", "veto_lncrna_smallrna", "veto_tad_risk_gene",
                  "veto_high_repeat_content", "pct_repeat",
                  "veto_ultraconserved_element", "max_50bp_rolling_phyloP",
                  "veto_external_regulatory_element", "veto_low_atac_accessibility", "no_rrbs_coverage",
                  "atac_mean", "atac_mean_percentile", "atac_mean_narrow", "atac_mean_narrow_percentile",
                  "atac_variability", "atac_peak_frequency",
                  "rrbs_mean", "rrbs_observed_fraction", "rrbs_observed_bp", "rrbs_variability",
                  "tad_boundary_distance", "mirna_distance", "risk_gene_radius_distance",
                  "gene_5prime_clearance", "gene_dense_clearance_old_bodydist",
                  "score_stability_atac", "score_low_peak_frequency", "score_low_methylation",
                  "score_stability_rrbs", "score_tad_distance", "final_score", "final_score_percentile"]

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
          f"{sum(1 for c in candidates if c['veto_tad_risk_gene'])} risk gene in own TAD, "
          f"{sum(1 for c in candidates if c['veto_high_repeat_content'])} high repeat content, "
          f"{sum(1 for c in candidates if c['veto_ultraconserved_element'])} ultraconserved element, "
          f"{sum(1 for c in candidates if c['veto_external_regulatory_element'])} external regulatory element overlap, "
          f"{sum(1 for c in candidates if c['veto_low_atac_accessibility'])} low ATAC accessibility "
          f"below p{args.min_atac_accessibility_percentile*100:.0f}), "
          f"{len(survivors)} ranked and passing.")
    print(f"Wrote full table to {args.out_scored}")
    print(f"Wrote ranked passing candidates BED to {args.out_passing_bed}")
    if survivors:
        top = survivors[0]
        print(f"Top candidate: {top['chrom']}:{top['start']}-{top['end']} "
              f"({top['orientation']}, score={top['final_score']}, "
              f"beats {top['final_score_percentile']}% of all 461 candidates)")


if __name__ == "__main__":
    main()
