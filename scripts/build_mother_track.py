#!/usr/bin/env python3
"""
Build the population "Mother Track" from the 71 per-dog CPM BigWigs and peak
sets, instead of intersecting 71 discrete peak calls (which zeroes out almost
everything). Produces three tracks:

  1. mother_track_weighted_mean.bw  - QC-weighted, local-confidence-weighted
     mean accessibility signal, binned genome-wide (the main Mother Track).
  2. peak_frequency.bw              - per merged-peak-region count of how many
     of the 71 dogs have their own peak there (context, not a filter).
  3. variability.bw                 - per-bin IQR of signal across the 71
     dogs (concordance metric, not a filter), computed only over dogs with
     meaningful local confidence at that bin.

A dog's raw (non-normalized) per-bin coverage is compared against its own
genome-wide background rate (lambda) to get a local "confidence" in [0,1]:
confidence = raw / (raw + lambda). This is what lets a low-depth dog's small
raw count still count as real signal, while the same raw count from a
high-depth dog (whose background lambda is much higher) is recognized as
noise-floor and down-weighted - and lets bins with essentially no coverage
for a given dog be excluded from that bin's average instead of dragging it
toward zero as if it were a real "closed chromatin" measurement.

Signal aggregation is done one chromosome at a time to bound memory.
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


def compute_lambda(raw_bw_path, total_genome_bp):
    """Genome-wide mean raw coverage for a dog: the expected per-bin
    background rate under a 'no real signal here' null model."""
    bw = pyBigWig.open(raw_bw_path)
    header = bw.header()
    bw.close()
    lam = header.get("sumData", 0.0) / total_genome_bp
    return max(lam, 1e-6)


def read_binned(bw, chrom, size, bin_size, n_bins):
    """Fast per-bin mean via a single bulk values() read + numpy reshape,
    instead of bw.stats(nBins=...) which is extremely slow for
    multi-million-bin queries (bin-by-bin overhead in the C extension)."""
    if chrom not in bw.chroms():
        return np.zeros(n_bins, dtype=np.float32)
    arr = bw.values(chrom, 0, size, numpy=True)
    arr = np.nan_to_num(arr, nan=0.0).astype(np.float32)
    pad = n_bins * bin_size - size
    if pad > 0:
        arr = np.pad(arr, (0, pad))
    return arr.reshape(n_bins, bin_size).mean(axis=1)


def build_weighted_mean_and_variability(samples, weights, per_dog_dir, chrom_sizes,
                                         bin_size, out_mean_bw, out_var_bw,
                                         confidence_floor=0.2):
    w = np.array(weights, dtype=np.float64)
    total_genome_bp = sum(size for _, size in chrom_sizes)

    mean_bw = pyBigWig.open(out_mean_bw, "w")
    var_bw = pyBigWig.open(out_var_bw, "w")
    mean_bw.addHeader(chrom_sizes)
    var_bw.addHeader(chrom_sizes)

    cpm_handles = {}
    raw_handles = {}
    lambdas = np.zeros(len(samples), dtype=np.float64)
    for i, s in enumerate(samples):
        cpm_handles[s] = pyBigWig.open(os.path.join(per_dog_dir, s, "bigwig", f"{s}.cpm.bw"))
        raw_path = os.path.join(per_dog_dir, s, "bigwig", f"{s}.raw.bw")
        raw_handles[s] = pyBigWig.open(raw_path)
        lambdas[i] = compute_lambda(raw_path, total_genome_bp)

    for chrom, size in chrom_sizes:
        n_bins = math.ceil(size / bin_size)
        cpm_mat = np.zeros((len(samples), n_bins), dtype=np.float32)
        raw_mat = np.zeros((len(samples), n_bins), dtype=np.float32)

        for i, s in enumerate(samples):
            cpm_mat[i, :] = read_binned(cpm_handles[s], chrom, size, bin_size, n_bins)
            raw_mat[i, :] = read_binned(raw_handles[s], chrom, size, bin_size, n_bins)

        # local confidence per dog per bin: raw / (raw + lambda_dog).
        # Low local coverage relative to that dog's own background -> low
        # confidence, regardless of whether "low" reads as 0 or a small
        # nonzero number - this is the piece a fixed raw-count threshold
        # can't capture, since the noise floor differs per dog by depth.
        confidence = raw_mat / (raw_mat + lambdas[:, None])

        w_local = w[:, None] * confidence
        denom = w_local.sum(axis=0)
        has_evidence = denom > 1e-9
        weighted_mean = np.zeros(n_bins, dtype=np.float64)
        weighted_mean[has_evidence] = (
            (cpm_mat * w_local).sum(axis=0)[has_evidence] / denom[has_evidence]
        )

        # variability: IQR computed only over dogs with real local evidence
        # (confidence above a floor) at that bin, so a bin isn't scored as
        # "highly variable" just because most dogs simply weren't measured
        # there.
        masked = np.where(confidence >= confidence_floor, cpm_mat, np.nan)
        with np.errstate(invalid="ignore", all="ignore"):
            q75 = np.nanpercentile(masked, 75, axis=0)
            q25 = np.nanpercentile(masked, 25, axis=0)
        variability = np.nan_to_num(q75 - q25, nan=0.0)

        starts = list(range(0, n_bins * bin_size, bin_size))[:n_bins]
        ends = [min(s + bin_size, size) for s in starts]
        chroms_col = [chrom] * n_bins

        mean_bw.addEntries(chroms_col, starts, ends=ends, values=[float(v) for v in weighted_mean])
        var_bw.addEntries(chroms_col, starts, ends=ends, values=[float(v) for v in variability])

        print(f"[weighted_mean/variability] {chrom}: {n_bins} bins done "
              f"(evidence in {int(has_evidence.sum())}/{n_bins} bins)")
        del cpm_mat, raw_mat, confidence, w_local
        gc.collect()

    for bw in cpm_handles.values():
        bw.close()
    for bw in raw_handles.values():
        bw.close()
    mean_bw.close()
    var_bw.close()


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
    ap.add_argument("--config", required=True, help="config.yaml with paths/params")
    ap.add_argument("--weights", required=True, help="qc_weights.tsv from compute_qc_weights.py")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    with open(args.config, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    work_dir = config["paths"]["work_dir"]
    per_dog_dir = os.path.join(work_dir, "per_dog")
    chrom_sizes_path = os.path.join(work_dir, "reference", "chrom.sizes")
    bin_size = config["params"]["bin_size"]

    samples, weights = load_weights(args.weights)
    chrom_sizes = load_chrom_sizes(chrom_sizes_path)

    os.makedirs(args.out_dir, exist_ok=True)
    mean_bw_path = os.path.join(args.out_dir, "mother_track_weighted_mean.bw")
    var_bw_path = os.path.join(args.out_dir, "variability.bw")
    freq_bedgraph_path = os.path.join(args.out_dir, "peak_frequency.bedgraph")
    freq_bw_path = os.path.join(args.out_dir, "peak_frequency.bw")

    print(f"Building weighted mean + variability across {len(samples)} dogs, bin={bin_size}bp")
    build_weighted_mean_and_variability(
        samples, weights, per_dog_dir, chrom_sizes, bin_size, mean_bw_path, var_bw_path
    )

    print("Building peak frequency track (bedtools multiinter)")
    build_peak_frequency(samples, per_dog_dir, chrom_sizes_path, freq_bedgraph_path, freq_bw_path)

    print("Done:")
    print(f"  {mean_bw_path}")
    print(f"  {var_bw_path}")
    print(f"  {freq_bw_path}")


if __name__ == "__main__":
    main()
