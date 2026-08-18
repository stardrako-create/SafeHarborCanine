#!/usr/bin/env python3
"""
Per-dog RRBS QC: mapping efficiency (Bismark PE report), an estimated
bisulfite non-conversion rate (from CHH/CHG methylation - real unmethylated
cytosines outside CpG context should convert almost completely, so residual
"methylation" there is mostly failed conversion, not biology), and CpG
coverage depth (from the .cov file) - the depth signal used later exactly
like ATAC's raw-coverage confidence weighting, since RRBS coverage is
extremely uneven (concentrated at MspI cut sites) and a fixed threshold
would misread "not covered here" as "unmethylated here".
"""
import argparse
import gzip
import re


def parse_report(path):
    """Generic 'Label: value' report parser - Bismark reports differ in
    exact wording across versions, so match by keyword rather than an exact
    line, and pull the first number found on a matching line."""
    values = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            key = line.lower()
            if "mapping efficiency" in key:
                target = "mapping_efficiency"
            elif "sequence pairs analysed" in key or "sequence pairs analyzed" in key:
                target = "pairs_analysed"
            elif "c methylated in cpg context" in key:
                target = "pct_meth_cpg"
            elif "c methylated in chg context" in key:
                target = "pct_meth_chg"
            elif "c methylated in chh context" in key:
                target = "pct_meth_chh"
            else:
                continue
            m = re.search(r"([\d.]+)\s*%?\s*$", line.strip())
            if not m:
                continue
            try:
                values[target] = float(m.group(1))
            except ValueError:
                continue
    return values


def coverage_stats(cov_path, min_coverage=5):
    total_positions = 0
    covered_positions = 0
    depth_sum = 0
    with gzip.open(cov_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            depth = int(cols[4]) + int(cols[5])
            total_positions += 1
            depth_sum += depth
            if depth >= min_coverage:
                covered_positions += 1
    mean_depth = depth_sum / total_positions if total_positions else 0.0
    return total_positions, covered_positions, mean_depth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--pe-report", required=True)
    ap.add_argument("--splitting-report", required=True)
    ap.add_argument("--cov", required=True)
    ap.add_argument("--min-coverage", type=int, default=5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pe_vals = parse_report(args.pe_report)
    split_vals = parse_report(args.splitting_report)

    mapping_eff = pe_vals.get("mapping_efficiency", 0.0)
    pct_chh = split_vals.get("pct_meth_chh", 0.0)
    pct_chg = split_vals.get("pct_meth_chg", 0.0)
    pct_cpg = split_vals.get("pct_meth_cpg", 0.0)
    non_conversion_proxy = (pct_chh + pct_chg) / 2.0
    conversion_rate = 100.0 - non_conversion_proxy

    n_cpg, n_covered, mean_depth = coverage_stats(args.cov, args.min_coverage)
    frac_covered = n_covered / n_cpg if n_cpg else 0.0

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("sample\tmapping_efficiency\tpct_meth_cpg\tpct_meth_chg\tpct_meth_chh\t"
                "bisulfite_conversion_rate\tn_cpg_called\tn_cpg_covered_min\t"
                "frac_cpg_covered\tmean_cpg_depth\n")
        f.write(f"{args.sample}\t{mapping_eff:.3f}\t{pct_cpg:.3f}\t{pct_chg:.3f}\t{pct_chh:.3f}\t"
                f"{conversion_rate:.3f}\t{n_cpg}\t{n_covered}\t{frac_covered:.4f}\t{mean_depth:.3f}\n")

    print(f"{args.sample}: mapping_eff={mapping_eff:.1f}% conversion~={conversion_rate:.1f}% "
          f"mean_CpG_depth={mean_depth:.2f}")


if __name__ == "__main__":
    main()
