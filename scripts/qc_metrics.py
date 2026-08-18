#!/usr/bin/env python3
"""
Per-sample QC metrics used later to weight each dog in the Mother Track:
  - FRiP: fraction of reads (in the final filtered BAM) that fall inside that
    dog's own MACS2 peaks.
  - TSS enrichment: mean read depth in +/-tss_flank windows around gene TSS,
    divided by mean read depth in the immediately flanking background windows.
    This is a simplified proxy for the ENCODE TSS enrichment score, not the
    exact ENCODE definition, but consistent across all 71 dogs which is what
    matters for relative weighting.
"""
import argparse
import os
import subprocess
import tempfile


def sh(cmd):
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def count_reads(bam, region_bed=None):
    cmd = ["samtools", "view", "-c"]
    if region_bed:
        cmd += ["-L", region_bed]
    cmd.append(bam)
    return int(sh(cmd).strip())


def bed_total_bp(bed_path):
    total = 0
    with open(bed_path, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            total += int(cols[2]) - int(cols[1])
    return total


def bedcov_sum(bed_path, bam_path):
    out = sh(["samtools", "bedcov", bed_path, bam_path])
    total = 0
    for line in out.splitlines():
        if not line.strip():
            continue
        total += int(line.split("\t")[-1])
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--peaks", required=True)
    ap.add_argument("--tss", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    ap.add_argument("--tss-flank", type=int, default=1000)
    ap.add_argument("--tss-bg", type=int, default=1000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    total_reads = count_reads(args.bam)
    reads_in_peaks = count_reads(args.bam, region_bed=args.peaks)
    frip = reads_in_peaks / total_reads if total_reads else 0.0

    with tempfile.TemporaryDirectory() as tmp:
        tss_window = os.path.join(tmp, "tss_window.bed")
        bg_window = os.path.join(tmp, "bg_window.bed")

        with open(tss_window, "w", encoding="utf-8") as f:
            f.write(sh(["bedtools", "slop", "-i", args.tss, "-g", args.chrom_sizes,
                        "-b", str(args.tss_flank)]))
        with open(bg_window, "w", encoding="utf-8") as f:
            f.write(sh(["bedtools", "flank", "-i", tss_window, "-g", args.chrom_sizes,
                        "-b", str(args.tss_bg)]))

        tss_bp = bed_total_bp(tss_window)
        bg_bp = bed_total_bp(bg_window)
        tss_depth_sum = bedcov_sum(tss_window, args.bam)
        bg_depth_sum = bedcov_sum(bg_window, args.bam)

    mean_tss_depth = tss_depth_sum / tss_bp if tss_bp else 0.0
    mean_bg_depth = bg_depth_sum / bg_bp if bg_bp else 0.0
    tss_enrichment = mean_tss_depth / (mean_bg_depth + 1e-9)

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("sample\ttotal_reads\treads_in_peaks\tfrip\tmean_tss_depth\tmean_bg_depth\ttss_enrichment\n")
        f.write(f"{args.sample}\t{total_reads}\t{reads_in_peaks}\t{frip:.6f}\t"
                f"{mean_tss_depth:.6f}\t{mean_bg_depth:.6f}\t{tss_enrichment:.6f}\n")

    print(f"{args.sample}: FRiP={frip:.4f} TSS_enrichment={tss_enrichment:.4f}")


if __name__ == "__main__":
    main()
