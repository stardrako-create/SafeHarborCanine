#!/usr/bin/env python3
"""
v1.20.0 unified Mother Track builder - one run across ALL dogs (any
cohort/batch), no post-hoc cohort joining. Replaces build_mother_track.py +
join_atac_cohorts.py for future runs (both kept, unused, for history).

Design, consolidated from what 2026-09-08/09's investigation established:

1. GATE, not soft-weight. A dog's per-bin confidence = raw/(raw+lambda)
   (lambda = that dog's own genome-wide mean raw coverage, its background
   rate). A dog below confidence_floor at a bin is excluded from that bin's
   mean entirely - not included at a softly-scaled-down weight. The earlier
   (pre-v1.20.0) version used confidence as a continuous multiplier with no
   hard cutoff, which let low-confidence dogs damp real peaks toward
   background - the Mother Track's dynamic range came out so compressed
   that no percentile cutoff on it could separate peaks from background at
   any threshold (0/32 confirmed real peaks reached even the 90th
   percentile). Above the floor, a dog counts at its own QC weight, no
   further continuous scaling by confidence magnitude.

2. GAIN: the gated mean is expressed as fold-enrichment over its own
   genome-wide background level (median across bins with evidence) -
   background lands near 1x, real signal shows as multiples above it - the
   units ATAC signal is normally discussed in, not a raw compressed CPM
   scale discovered by accident.

3. This track (renamed here: mother_track_accessibility_level.bw) is NOT
   the primary peak-discriminative signal - empirically verified
   2026-09-08 that even gated+gain-corrected, confirmed real peaks land at
   only the ~83rd percentile here. It answers "how open is chromatin",
   used only for a coarse accessibility floor. peak_frequency.bw (vote
   count: how many dogs have their own independently-called MACS3 peak
   nearby) is the actual discriminative track, unaffected by any of this -
   built the same way as before, see build_peak_frequency() below.

4. No raw.bw needed. lambda and total_reads are read from each dog's
   qc.tsv (persisted there by qc_metrics.py / backfill_lambda_v2.py); a
   dog's per-bin raw coverage is reconstructed from its cpm.bw as
   raw(bin) = cpm(bin) * total_reads / 1e6 - verified exact (ratio 1.0000
   at several spot-checked bins, 2026-09-09) against the actual raw.bw
   values before this script was trusted to drop raw.bw as an input.

5. Fixed 2026-09-10: variability.bw no longer writes 0.0 for bins where
   0 or 1 dogs passed the confidence gate. An IQR needs at least 2 values
   to mean anything - with 0, np.nanpercentile over an all-NaN slice was
   silently nan_to_num'd to 0.0; with 1, the IQR of a single value really
   is 0, but is equally uninformative about stability across dogs. Both
   were indistinguishable downstream from genuine stability (the highest
   possible score_stability_atac), the opposite of what they actually
   mean: no/insufficient evidence. Now: bins with fewer than
   min_variability_evidence_dogs (default 2) passing the gate are left as
   a real gap in variability.bw, matching the sparse-means-missing
   convention already used for peak_frequency.bw - see
   bw_utils.keep_bins_with_evidence().
"""
import argparse
import csv
import gc
import math
import os
import subprocess
import tempfile

import numpy as np
import pyBigWig
import yaml

from bw_utils import keep_bins_with_evidence


def load_chrom_sizes(path):
    chroms = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            name, size = line.rstrip("\n").split("\t")
            chroms.append((name, int(size)))
    return chroms


def load_weights(path):
    samples, weights = [], []
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            samples.append(row["sample"])
            weights.append(float(row["weight"]))
    return samples, weights


def load_lambda_and_scale(per_dog_dir, samples):
    """lambda and the CPM scale factor (total_reads/1e6) per dog, from
    qc.tsv - no raw.bw opened."""
    lambdas, scales = [], []
    for s in samples:
        qc_path = os.path.join(per_dog_dir, s, "qc", f"{s}.qc.tsv")
        with open(qc_path, encoding="utf-8") as f:
            row = next(csv.DictReader(f, delimiter="\t"))
        if "lambda" not in row:
            raise SystemExit(
                f"{s}: qc.tsv has no lambda column - run backfill_lambda_v2.py first"
            )
        lambdas.append(float(row["lambda"]))
        scales.append(float(row["total_reads"]) / 1e6)
    return np.array(lambdas, dtype=np.float64), np.array(scales, dtype=np.float64)


def read_binned(bw, chrom, size, bin_size, n_bins, sample_name="?", chunk_size=1_000_000):
    if chrom not in bw.chroms():
        return np.zeros(n_bins, dtype=np.float32)
    try:
        arr = bw.values(chrom, 0, size, numpy=True)
    except RuntimeError:
        # Found 2026-09-11: a "whole chromosome unreadable" error can come
        # from ONE bad data block on an otherwise-fine chromosome, not
        # genuine corruption everywhere - confirmed directly on
        # SRR27376673/NC_051811.1 (RIT2's own window there read perfectly
        # fine once queried in isolation; only a different sub-range
        # elsewhere on the chromosome was bad). Falling back to zero for
        # the WHOLE chromosome on any failure (as the previous version of
        # this fix did) throws away real, readable data unnecessarily.
        # Localize instead: retry chunk-by-chunk, zero only the specific
        # chunks that individually fail.
        arr = np.zeros(size, dtype=np.float32)
        for chunk_start in range(0, size, chunk_size):
            chunk_end = min(chunk_start + chunk_size, size)
            try:
                arr[chunk_start:chunk_end] = bw.values(chrom, chunk_start, chunk_end, numpy=True)
            except RuntimeError:
                print(f"  [warn] {sample_name}: could not read {chrom}:{chunk_start}-{chunk_end} "
                      f"- zero-contribution for this ~{chunk_size:,}bp chunk only")
    arr = np.nan_to_num(arr, nan=0.0).astype(np.float32)
    pad = n_bins * bin_size - size
    if pad > 0:
        arr = np.pad(arr, (0, pad))
    return arr.reshape(n_bins, bin_size).mean(axis=1)


def build_accessibility_and_variability(samples, weights, lambdas, scales, per_dog_dir,
                                         chrom_sizes, bin_size, out_mean_bw, out_var_bw,
                                         confidence_floor=0.2, min_variability_evidence_dogs=2):
    w = np.array(weights, dtype=np.float64)

    var_bw = pyBigWig.open(out_var_bw, "w")
    var_bw.addHeader(chrom_sizes)

    cpm_handles = {s: pyBigWig.open(os.path.join(per_dog_dir, s, "bigwig", f"{s}.cpm.bw"))
                   for s in samples}

    buffered_mean, buffered_evidence = {}, {}
    for chrom, size in chrom_sizes:
        n_bins = math.ceil(size / bin_size)
        cpm_mat = np.zeros((len(samples), n_bins), dtype=np.float32)
        for i, s in enumerate(samples):
            try:
                cpm_mat[i, :] = read_binned(cpm_handles[s], chrom, size, bin_size, n_bins,
                                             sample_name=s)
            except RuntimeError as e:
                # read_binned() already retries chunk-by-chunk and zeros only
                # the specific chunks that fail (found 2026-09-11: a "whole
                # chromosome unreadable" error can come from one bad block,
                # not genuine corruption everywhere - confirmed directly on
                # SRR27376673/NC_051811.1, where RIT2's own window read fine
                # in isolation). This outer catch is now just a defensive
                # fallback for something read_binned() itself can't recover
                # from; cpm_mat's row is already zero-initialized.
                print(f"  [warn] {s}: could not read {chrom} at all ({e}) - "
                      f"treating as zero-contribution for {chrom} only")

        raw_mat = cpm_mat * scales[:, None]
        confidence = raw_mat / (raw_mat + lambdas[:, None])

        gate = confidence >= confidence_floor
        w_gated = np.where(gate, w[:, None], 0.0)
        denom = w_gated.sum(axis=0)
        has_evidence = denom > 1e-9
        weighted_mean = np.zeros(n_bins, dtype=np.float64)
        weighted_mean[has_evidence] = (
            (cpm_mat * w_gated).sum(axis=0)[has_evidence] / denom[has_evidence]
        )

        masked = np.where(gate, cpm_mat, np.nan)
        with np.errstate(invalid="ignore", all="ignore"):
            q75 = np.nanpercentile(masked, 75, axis=0)
            q25 = np.nanpercentile(masked, 25, axis=0)
        variability = np.nan_to_num(q75 - q25, nan=0.0)
        n_gated = gate.sum(axis=0)
        has_var_evidence = keep_bins_with_evidence(n_gated, min_variability_evidence_dogs)

        starts = np.arange(0, n_bins * bin_size, bin_size, dtype=np.int64)[:n_bins]
        ends = np.minimum(starts + bin_size, size)
        n_var_entries = int(has_var_evidence.sum())
        if n_var_entries > 0:
            # pyBigWig's addEntries rejects a genuinely empty call outright
            # (confirmed 2026-09-10: raises "out of order, precede already
            # added entries" even in isolation) - small/sparse scaffolds can
            # legitimately have zero bins clearing the evidence threshold.
            var_bw.addEntries(
                [chrom] * n_var_entries,
                starts[has_var_evidence].tolist(),
                ends=ends[has_var_evidence].tolist(),
                values=[float(v) for v in variability[has_var_evidence]],
            )

        buffered_mean[chrom] = weighted_mean.astype(np.float32)
        buffered_evidence[chrom] = has_evidence

        print(f"[accessibility/variability] {chrom}: {n_bins} bins done "
              f"(evidence in {int(has_evidence.sum())}/{n_bins} bins; "
              f"variability written for {int(has_var_evidence.sum())}/{n_bins} bins, "
              f"{n_bins - int(has_var_evidence.sum())} skipped for "
              f"<{min_variability_evidence_dogs} dogs)")
        del cpm_mat, raw_mat, confidence, gate, w_gated, masked, weighted_mean, variability, n_gated, has_var_evidence
        gc.collect()

    for bw in cpm_handles.values():
        bw.close()
    var_bw.close()

    all_means = np.concatenate([buffered_mean[c] for c, _ in chrom_sizes])
    all_evidence = np.concatenate([buffered_evidence[c] for c, _ in chrom_sizes])
    background = float(np.median(all_means[all_evidence])) if all_evidence.any() else 1.0
    background = max(background, 1e-9)
    evidence_vals = all_means[all_evidence]
    print(f"[gain] genome-wide background level (median of gated weighted mean): {background:.6g}")
    print(
        "[gain] pre-rescale gated weighted-mean percentiles: "
        f"p50={np.percentile(evidence_vals, 50):.4g} "
        f"p90={np.percentile(evidence_vals, 90):.4g} "
        f"p99={np.percentile(evidence_vals, 99):.4g} "
        f"p99.9={np.percentile(evidence_vals, 99.9):.4g} "
        f"max={evidence_vals.max():.4g}"
    )
    del all_means, all_evidence, evidence_vals
    gc.collect()

    mean_bw = pyBigWig.open(out_mean_bw, "w")
    mean_bw.addHeader(chrom_sizes)
    for chrom, size in chrom_sizes:
        n_bins = math.ceil(size / bin_size)
        fold = buffered_mean[chrom] / background
        starts = list(range(0, n_bins * bin_size, bin_size))[:n_bins]
        ends = [min(st + bin_size, size) for st in starts]
        chroms_col = [chrom] * n_bins
        mean_bw.addEntries(chroms_col, starts, ends=ends, values=[float(v) for v in fold])
    mean_bw.close()


def sorted_bed_from_narrowpeak(narrowpeak_path, out_path):
    intervals = []
    with open(narrowpeak_path, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            intervals.append((cols[0], int(cols[1]), int(cols[2])))
    intervals.sort(key=lambda t: (t[0], t[1]))
    with open(out_path, "w", encoding="utf-8") as f:
        for chrom, s, e in intervals:
            f.write(f"{chrom}\t{s}\t{e}\n")


def build_peak_frequency(samples, per_dog_dir, chrom_sizes_path, out_bedgraph, out_bw):
    with tempfile.TemporaryDirectory() as tmp:
        sorted_beds = []
        for s in samples:
            narrowpeak = os.path.join(per_dog_dir, s, "peaks", f"{s}_peaks.narrowPeak")
            sorted_path = os.path.join(tmp, f"{s}.sorted.bed")
            sorted_bed_from_narrowpeak(narrowpeak, sorted_path)
            sorted_beds.append(sorted_path)

        result = subprocess.run(
            ["bedtools", "multiinter", "-i"] + sorted_beds,
            check=True, capture_output=True, text=True,
        )

    with open(out_bedgraph, "w", encoding="utf-8") as f:
        for line in result.stdout.splitlines():
            cols = line.split("\t")
            chrom, start, end, num_samples = cols[0], cols[1], cols[2], cols[3]
            f.write(f"{chrom}\t{start}\t{end}\t{num_samples}\n")

    subprocess.run(["bedGraphToBigWig", out_bedgraph, chrom_sizes_path, out_bw], check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True, nargs="+", help="unified sample list, any cohort")
    ap.add_argument("--weights", required=True, help="compute_qc_weights.py output, unified across all samples")
    ap.add_argument("--per-dog-dir", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--bin-size", type=int, default=25)
    ap.add_argument("--confidence-floor", type=float, default=0.2)
    ap.add_argument("--min-variability-evidence-dogs", type=int, default=2,
                     help="bins where fewer dogs pass the confidence gate are left as "
                          "a gap in variability.bw instead of a false 0.0 'stable'")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    weight_samples, weights = load_weights(args.weights)
    if set(weight_samples) != set(args.samples):
        raise SystemExit(
            f"--samples ({len(args.samples)}) and --weights ({len(weight_samples)}) "
            "sample sets don't match - weights must be computed across exactly the "
            "unified sample list this run uses"
        )
    # order samples exactly as the weights file lists them, so weights[i] matches samples[i]
    samples = weight_samples

    lambdas, scales = load_lambda_and_scale(args.per_dog_dir, samples)
    chrom_sizes = load_chrom_sizes(args.chrom_sizes)

    os.makedirs(args.out_dir, exist_ok=True)
    mean_bw_path = os.path.join(args.out_dir, "mother_track_accessibility_level.bw")
    var_bw_path = os.path.join(args.out_dir, "variability.bw")
    freq_bedgraph_path = os.path.join(args.out_dir, "peak_frequency.bedgraph")
    freq_bw_path = os.path.join(args.out_dir, "peak_frequency.bw")

    print(f"Building accessibility-level + variability across {len(samples)} dogs "
          f"(unified, no cohort split), bin={args.bin_size}bp")
    build_accessibility_and_variability(
        samples, weights, lambdas, scales, args.per_dog_dir, chrom_sizes,
        args.bin_size, mean_bw_path, var_bw_path,
        confidence_floor=args.confidence_floor,
        min_variability_evidence_dogs=args.min_variability_evidence_dogs,
    )

    print("Building peak frequency track (bedtools multiinter)")
    build_peak_frequency(samples, args.per_dog_dir, args.chrom_sizes, freq_bedgraph_path, freq_bw_path)

    print("Done:")
    print(f"  {mean_bw_path}")
    print(f"  {var_bw_path}")
    print(f"  {freq_bw_path}")


if __name__ == "__main__":
    main()
