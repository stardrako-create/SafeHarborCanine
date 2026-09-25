#!/usr/bin/env python3
"""
Combine two independently-built ATAC Mother Tracks (71-dog + 5-dog Ehsan
cohorts) into one joined track, by cohort size rather than re-pooling all
76 raw dogs as equally-weighted individuals:

    joined(bin) = [ n_a * mean_a(bin) + n_b * mean_b(bin) ] / (n_a + n_b)

Applied per-bin to both the weighted-mean and variability tracks (see
04_tracks_processadas/ROS_Cfam_1.0/METHODS.md, "Joined 76-dog track" - same
formula, re-run here against the gate+gain-fixed cohort tracks).
"""
import argparse

import pyBigWig


def join_bigwigs(bw_a_path, n_a, bw_b_path, n_b, out_path):
    bw_a = pyBigWig.open(bw_a_path)
    bw_b = pyBigWig.open(bw_b_path)
    chrom_sizes = list(bw_a.chroms().items())
    total_n = n_a + n_b

    out_bw = pyBigWig.open(out_path, "w")
    out_bw.addHeader(chrom_sizes)

    for chrom, size in chrom_sizes:
        ivs_a = bw_a.intervals(chrom) or []
        ivs_b = bw_b.intervals(chrom) or []
        if len(ivs_a) != len(ivs_b):
            raise ValueError(
                f"{chrom}: bin count mismatch between cohorts "
                f"({len(ivs_a)} vs {len(ivs_b)}) - both tracks must share bin_size"
            )
        starts = [s for s, _, _ in ivs_a]
        ends = [e for _, e, _ in ivs_a]
        joined = [
            (n_a * va + n_b * vb) / total_n
            for (_, _, va), (_, _, vb) in zip(ivs_a, ivs_b)
        ]
        out_bw.addEntries([chrom] * len(starts), starts, ends=ends, values=joined)
        print(f"[join] {chrom}: {len(starts)} bins joined")

    bw_a.close()
    bw_b.close()
    out_bw.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mean-a", required=True)
    ap.add_argument("--mean-b", required=True)
    ap.add_argument("--var-a", required=True)
    ap.add_argument("--var-b", required=True)
    ap.add_argument("--n-a", type=int, required=True)
    ap.add_argument("--n-b", type=int, required=True)
    ap.add_argument("--out-mean", required=True)
    ap.add_argument("--out-var", required=True)
    args = ap.parse_args()

    print(f"Joining weighted mean ({args.n_a} + {args.n_b} dogs)...")
    join_bigwigs(args.mean_a, args.n_a, args.mean_b, args.n_b, args.out_mean)
    print(f"Joining variability ({args.n_a} + {args.n_b} dogs)...")
    join_bigwigs(args.var_a, args.n_a, args.var_b, args.n_b, args.out_var)
    print("Done:")
    print(f"  {args.out_mean}")
    print(f"  {args.out_var}")


if __name__ == "__main__":
    main()
