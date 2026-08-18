#!/usr/bin/env python3
"""
Build the population "Mother Track" for RRBS methylation from the 71 per-dog
Bismark .cov.gz files, analogous to build_mother_track.py for ATAC. RRBS
coverage is a sparse list of CpG positions (reduced representation,
concentrated near MspI cut sites) rather than a continuous per-base bigwig,
so each dog's .cov.gz is parsed once (vectorized) into whole-genome bins and
folded into running totals - one dog at a time, so peak memory never holds
more than ~2 dogs' worth of per-bin arrays plus the running totals.

Two weight layers, same philosophy as the ATAC track:
  - GLOBAL per-dog weight w_i: from bisulfite_conversion_rate and
    mapping_efficiency in {sample}.rrbs_qc.tsv (written by qc_metrics_rrbs.py),
    min-max normalized and averaged, rescaled onto [floor, 1.0] so no dog is
    ever fully zeroed out - same recipe as compute_qc_weights.py for ATAC.
  - LOCAL per-bin confidence per dog: bin_depth / (bin_depth + lambda_i),
    where lambda_i is that dog's own mean_cpg_depth (already computed in its
    qc.tsv) - the RRBS analogue of ATAC's raw/(raw+lambda) trick. It's needed
    because RRBS depth is extremely uneven (concentrated at MspI sites), so a
    fixed coverage threshold would misread "not covered here" as "0%
    methylated here" instead of "not measured here".

Produces three genome-wide BigWigs (all dense - every bin gets a value, 0
where there is no evidence, matching build_mother_track.py's convention):
  1. methylation_weighted_mean.bw - QC-weighted, local-confidence-weighted
     mean %CpG methylation per bin (the main Mother Track), expressed as a
     percentage (0-100), matching the pct_meth_cpg convention already used
     in {sample}.rrbs_qc.tsv.
  2. cpg_coverage_frequency.bw   - per-bin count of dogs with meaningful
     local confidence there (context, not a filter).
  3. variability.bw              - per-bin standard deviation of %methylation
     across dogs with real local evidence (concordance metric, not a
     filter). Standard deviation (not IQR, unlike the ATAC script) because it
     can be accumulated with running sums in the same single per-dog pass
     used for the mean track, instead of requiring all 71 dogs' per-bin
     values held in memory simultaneously across the whole genome.
"""
import argparse
import csv
import gc
import math
import os

import numpy as np
import pandas as pd
import pyBigWig
import yaml


def load_chrom_sizes(path):
    chroms = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            name, size = line.rstrip("\n").split("\t")
            chroms.append((name, int(size)))
    return chroms


def minmax(values):
    lo, hi = min(values), max(values)
    if hi - lo < 1e-12:
        return [1.0 for _ in values]
    return [(v - lo) / (hi - lo) for v in values]


def load_qc(per_dog_dir, samples):
    """Read each sample's {sample}.rrbs_qc.tsv (one data row) for the global
    weight inputs and the per-dog lambda (mean_cpg_depth)."""
    rows = {}
    for s in samples:
        path = os.path.join(per_dog_dir, s, "qc", f"{s}.rrbs_qc.tsv")
        with open(path, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows[s] = next(reader)
    return rows


def compute_global_weights(qc_rows, samples, floor):
    mapping_eff = [float(qc_rows[s]["mapping_efficiency"]) for s in samples]
    conversion = [float(qc_rows[s]["bisulfite_conversion_rate"]) for s in samples]
    me_norm = minmax(mapping_eff)
    cv_norm = minmax(conversion)
    weights = []
    for mn, cn in zip(me_norm, cv_norm):
        composite = (mn + cn) / 2.0
        weights.append(floor + (1.0 - floor) * composite)
    return weights


def bin_one_dog(cov_path, chrom_names, chrom_dtype, offsets, nbins_arr, bin_size, total_bins):
    """Vectorized parse + genome-wide binning of one dog's .cov.gz. Returns
    (depth_binned, meth_binned), each a float64 array of length total_bins."""
    df = pd.read_csv(
        cov_path, sep="\t", header=None,
        names=["chrom", "pos", "end", "pctmeth", "meth", "unmeth"],
        usecols=["chrom", "pos", "meth", "unmeth"],
        dtype={"pos": "int64", "meth": "int32", "unmeth": "int32"},
        compression="gzip",
    )
    df["chrom"] = df["chrom"].astype(chrom_dtype)
    codes = df["chrom"].cat.codes.to_numpy()
    valid = codes >= 0
    if not valid.all():
        dropped = int((~valid).sum())
        print(f"  [warn] {dropped} rows on chromosomes outside chrom.sizes, skipping")
    codes = codes[valid]
    pos = df["pos"].to_numpy()[valid]
    meth = df["meth"].to_numpy()[valid]
    unmeth = df["unmeth"].to_numpy()[valid]

    local_bin = np.minimum((pos - 1) // bin_size, nbins_arr[codes] - 1)
    global_bin = offsets[codes] + local_bin

    depth = (meth + unmeth).astype(np.float64)
    depth_binned = np.bincount(global_bin, weights=depth, minlength=total_bins)
    meth_binned = np.bincount(global_bin, weights=meth.astype(np.float64), minlength=total_bins)
    return depth_binned, meth_binned


def build_tracks(samples, weights, qc_rows, per_dog_dir, chrom_sizes, bin_size,
                  out_mean_bw, out_freq_bw, out_var_bw, confidence_floor):
    chrom_names = [c for c, _ in chrom_sizes]
    chrom_dtype = pd.CategoricalDtype(categories=chrom_names, ordered=True)
    nbins_per_chrom = [math.ceil(size / bin_size) for _, size in chrom_sizes]
    nbins_arr = np.array(nbins_per_chrom, dtype=np.int64)
    offsets = np.concatenate(([0], np.cumsum(nbins_arr)))[:-1]
    total_bins = int(nbins_arr.sum())
    print(f"Genome-wide bins: {total_bins:,} (bin_size={bin_size}bp, "
          f"{len(chrom_sizes)} contigs)")

    mean_num = np.zeros(total_bins, dtype=np.float64)
    mean_den = np.zeros(total_bins, dtype=np.float64)
    freq_count = np.zeros(total_bins, dtype=np.float64)
    var_sum = np.zeros(total_bins, dtype=np.float64)
    var_sumsq = np.zeros(total_bins, dtype=np.float64)
    var_n = np.zeros(total_bins, dtype=np.float64)

    for i, s in enumerate(samples):
        lam = max(float(qc_rows[s]["mean_cpg_depth"]), 1e-3)
        cov_path = os.path.join(per_dog_dir, s, "meth", f"{s}_pe.bismark.cov.gz")
        depth_binned, meth_binned = bin_one_dog(
            cov_path, chrom_names, chrom_dtype, offsets, nbins_arr, bin_size, total_bins
        )

        confidence = depth_binned / (depth_binned + lam)
        pct_meth = np.divide(
            meth_binned, depth_binned,
            out=np.zeros_like(meth_binned), where=depth_binned > 0,
        )

        w_i = weights[i]
        mean_num += w_i * confidence * pct_meth
        mean_den += w_i * confidence

        evidence = confidence >= confidence_floor
        freq_count += evidence
        var_sum += np.where(evidence, pct_meth, 0.0)
        var_sumsq += np.where(evidence, pct_meth * pct_meth, 0.0)
        var_n += evidence

        print(f"[{i + 1}/{len(samples)}] {s}: weight={w_i:.3f} lambda={lam:.2f} "
              f"bins_with_depth={int((depth_binned > 0).sum()):,}")

        del depth_binned, meth_binned, confidence, pct_meth, evidence
        gc.collect()

    has_evidence = mean_den > 1e-9
    weighted_mean = np.zeros(total_bins, dtype=np.float64)
    weighted_mean[has_evidence] = mean_num[has_evidence] / mean_den[has_evidence]
    weighted_mean *= 100.0  # fraction -> percentage, matches pct_meth_cpg convention

    has_var = var_n > 1
    var_mean = np.zeros(total_bins, dtype=np.float64)
    var_mean[var_n > 0] = var_sum[var_n > 0] / var_n[var_n > 0]
    variance = np.zeros(total_bins, dtype=np.float64)
    variance[has_var] = np.clip(
        var_sumsq[has_var] / var_n[has_var] - var_mean[has_var] ** 2, 0, None
    )
    variability = np.sqrt(variance) * 100.0  # same units as the mean track

    print("Writing BigWigs...")
    mean_bw = pyBigWig.open(out_mean_bw, "w")
    freq_bw = pyBigWig.open(out_freq_bw, "w")
    var_bw = pyBigWig.open(out_var_bw, "w")
    mean_bw.addHeader(chrom_sizes)
    freq_bw.addHeader(chrom_sizes)
    var_bw.addHeader(chrom_sizes)

    for ci, (chrom, size) in enumerate(chrom_sizes):
        n_bins = nbins_per_chrom[ci]
        lo = offsets[ci]
        hi = lo + n_bins
        starts = list(range(0, n_bins * bin_size, bin_size))[:n_bins]
        ends = [min(st + bin_size, size) for st in starts]
        chroms_col = [chrom] * n_bins

        mean_bw.addEntries(chroms_col, starts, ends=ends,
                            values=[float(v) for v in weighted_mean[lo:hi]])
        freq_bw.addEntries(chroms_col, starts, ends=ends,
                            values=[float(v) for v in freq_count[lo:hi]])
        var_bw.addEntries(chroms_col, starts, ends=ends,
                           values=[float(v) for v in variability[lo:hi]])

    mean_bw.close()
    freq_bw.close()
    var_bw.close()
    print(f"  evidence in {int(has_evidence.sum()):,}/{total_bins:,} bins genome-wide")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="config_rrbs.yaml with paths/params")
    ap.add_argument("--out-dir", default=None,
                     help="defaults to paths.final_dir in the config")
    ap.add_argument("--bin-size", type=int, default=50)
    ap.add_argument("--min-weight-floor", type=float, default=0.2)
    ap.add_argument("--confidence-floor", type=float, default=0.2)
    args = ap.parse_args()

    with open(args.config, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    samples = config["samples"]
    work_dir = config["paths"]["work_dir"]
    per_dog_dir = os.path.join(work_dir, "per_dog")
    ref_dir = config["paths"]["ref_dir"]
    chrom_sizes_path = os.path.join(ref_dir, "chrom.sizes")
    out_dir = args.out_dir or config["paths"]["final_dir"]

    print(f"Loading QC for {len(samples)} dogs...")
    qc_rows = load_qc(per_dog_dir, samples)
    weights = compute_global_weights(qc_rows, samples, args.min_weight_floor)
    chrom_sizes = load_chrom_sizes(chrom_sizes_path)

    os.makedirs(out_dir, exist_ok=True)
    weights_path = os.path.join(out_dir, "rrbs_qc_weights.tsv")
    with open(weights_path, "w", encoding="utf-8") as f:
        f.write("sample\tmapping_efficiency\tbisulfite_conversion_rate\tmean_cpg_depth\tweight\n")
        for s, w in zip(samples, weights):
            r = qc_rows[s]
            f.write(f"{s}\t{r['mapping_efficiency']}\t{r['bisulfite_conversion_rate']}\t"
                    f"{r['mean_cpg_depth']}\t{w:.6f}\n")
    print(f"Wrote weights for {len(samples)} dogs to {weights_path}")
    print(f"weight range: {min(weights):.3f} - {max(weights):.3f}")

    mean_bw_path = os.path.join(out_dir, "methylation_weighted_mean.bw")
    freq_bw_path = os.path.join(out_dir, "cpg_coverage_frequency.bw")
    var_bw_path = os.path.join(out_dir, "variability.bw")

    print(f"Building weighted mean + coverage frequency + variability across "
          f"{len(samples)} dogs, bin={args.bin_size}bp")
    build_tracks(samples, weights, qc_rows, per_dog_dir, chrom_sizes, args.bin_size,
                 mean_bw_path, freq_bw_path, var_bw_path,
                 confidence_floor=args.confidence_floor)

    print("Done:")
    print(f"  {mean_bw_path}")
    print(f"  {freq_bw_path}")
    print(f"  {var_bw_path}")
    print(f"  {weights_path}")


if __name__ == "__main__":
    main()
