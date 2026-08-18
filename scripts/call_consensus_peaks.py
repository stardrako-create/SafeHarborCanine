#!/usr/bin/env python3
"""
Phase 3: call consensus peaks directly from the peak_frequency track (how
many of the 71 dogs have their own individually-called MACS3 peak at a given
region), following the same logic ArchR's addReproduciblePeakSet() uses for
population-scale ATAC-seq: a peak is "reproducible" if it is present in a
data-size-dependent fraction of samples (ArchR's own documented example is
reproducibility="(n+1)/2", i.e. a simple majority). A percentile cutoff on
the continuous weighted-mean signal was tried first and rejected: that
track's nonzero-value distribution is narrow and non-bimodal (p50=0.09,
p99.9=0.83), so a magnitude threshold has no natural peak/background
separation to exploit - vote count across dogs does.

Each resulting consensus region is then annotated with the mean weighted
signal and mean variability from the Mother Track, as context - not as a
filter, matching the project's stated design (EpiC Dog/TAD-style annotation,
never an isolated exclusion criterion).
"""
import argparse
import os
import subprocess

import pyBigWig


def filter_and_merge(freq_bedgraph_path, min_count, max_gap, out_path):
    filtered_path = out_path + ".filtered.tmp"
    n_in = 0
    with open(freq_bedgraph_path, encoding="utf-8") as fin, \
         open(filtered_path, "w", encoding="utf-8") as fout:
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            if int(cols[3]) >= min_count:
                fout.write(line)
                n_in += 1

    subprocess.run(
        ["bedtools", "merge", "-d", str(max_gap), "-c", "4,4", "-o", "max,mean",
         "-i", filtered_path],
        check=True, stdout=open(out_path, "w", encoding="utf-8"),
    )
    os.remove(filtered_path)
    return n_in


def annotate_peaks(merged_path, mean_bw_path, var_bw_path, out_path, n_samples):
    mean_bw = pyBigWig.open(mean_bw_path)
    var_bw = pyBigWig.open(var_bw_path)

    n = 0
    with open(merged_path, encoding="utf-8") as fin, \
         open(out_path, "w", encoding="utf-8") as fout:
        fout.write("chrom\tstart\tend\tname\tmax_dogs_supporting\tmean_dogs_supporting\t"
                    "frac_dogs_supporting\tmean_weighted_signal\tmean_variability\n")
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            max_count, mean_count = float(cols[3]), float(cols[4])
            n += 1
            name = f"consensus_peak_{n}"

            sig_val = 0.0
            if chrom in mean_bw.chroms():
                stats = mean_bw.stats(chrom, start, end, type="mean")
                sig_val = stats[0] if stats and stats[0] is not None else 0.0

            var_val = 0.0
            if chrom in var_bw.chroms():
                stats = var_bw.stats(chrom, start, end, type="mean")
                var_val = stats[0] if stats and stats[0] is not None else 0.0

            frac = mean_count / n_samples
            fout.write(f"{chrom}\t{start}\t{end}\t{name}\t{max_count:.0f}\t{mean_count:.2f}\t"
                       f"{frac:.3f}\t{sig_val:.6f}\t{var_val:.6f}\n")

    mean_bw.close()
    var_bw.close()
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq-bedgraph", required=True)
    ap.add_argument("--mean-bw", required=True)
    ap.add_argument("--var-bw", required=True)
    ap.add_argument("--n-samples", type=int, required=True)
    ap.add_argument("--min-count", type=int, default=None,
                     help="Minimum number of dogs supporting a region. "
                          "Default: majority, ceil((n+1)/2), matching ArchR's "
                          "documented reproducibility=\"(n+1)/2\" convention.")
    ap.add_argument("--max-gap", type=int, default=75)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    min_count = args.min_count
    if min_count is None:
        min_count = -(-(args.n_samples + 1) // 2)  # ceil((n+1)/2)

    os.makedirs(args.out_dir, exist_ok=True)
    merged_path = os.path.join(args.out_dir, "consensus_peaks_merged.bed")
    final_path = os.path.join(args.out_dir, "consensus_peaks_ATAC.bed")

    print(f"Consensus threshold: >= {min_count}/{args.n_samples} dogs "
          f"({100*min_count/args.n_samples:.1f}%), max_gap={args.max_gap}bp")

    n_in = filter_and_merge(args.freq_bedgraph, min_count, args.max_gap, merged_path)
    print(f"{n_in} raw intervals passed the threshold before merging")

    n_out = annotate_peaks(merged_path, args.mean_bw, args.var_bw, final_path, args.n_samples)
    print(f"Done: {n_out} consensus peaks -> {final_path}")


if __name__ == "__main__":
    main()
