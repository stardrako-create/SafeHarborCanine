#!/usr/bin/env python3
"""
Compute CpG islands natively on ROS_Cfam_1.0 from the reference FASTA,
written 2026-09-11 to replace dependence on an opaque third-party
regulatory BED.

Why this exists: the project's only regulatory-element annotation was
05_SHIP/ehsan_regulatory_elements_ROS.bed - 75,600 bare intervals, no
element type, no evidence, no name, no header, received 2026-08-21 as a
CanFam3.1 set lifted to ROS_Cfam_1.0 by Ehsan. It cannot be audited (the
metadata was never in it), it carries whatever error that liftover
introduced, and it is already known to be out of date relative to Ehsan's
own current curation (his 2026-09-11 email cites a regulatory element at
NC_051807.1:10,779,807-10,837,318 that this file does not contain).

CpG islands were part of Ehsan's own methodology too ("CanFam3.1
regulatory elements/CpG islands liftover"), so this is comparable to his
approach rather than orthogonal to it - but computed directly on the
target assembly, so there is no liftover step to go wrong, and every
criterion is explicit and reproducible here rather than inherited.

Emits BOTH standard parameterizations, labelled, with the statistics that
justify each call:
  - CGI_GGF  Gardiner-Garden & Frommer 1987: GC >= 50%, Obs/Exp >= 0.6,
             length >= 200bp. The classic, permissive definition.
  - CGI_TJ   Takai & Jones 2002: GC >= 55%, Obs/Exp >= 0.65, length >=
             500bp. Stricter, designed specifically to exclude Alu and
             other repetitive false positives.
Each definition is emitted independently with its own coordinates and statistics.
Overlapping GGF and TJ calls are retained; no overlap-derived dual labels.

Obs/Exp = (n_CpG * window_length) / (n_C * n_G), the standard formula.
Windows below --min-acgt-frac non-N bases are skipped - assembly gaps
would otherwise produce meaningless ratios.

Output BED (tab-separated, with a # header line so it is self-describing,
unlike the file it replaces):
  chrom  start  end  type  length  gc_frac  obs_exp  n_cpg
"""
import argparse

import numpy as np
import pysam

PARAMS = {
    "CGI_GGF": dict(min_gc=0.50, min_oe=0.60, min_len=200),
    "CGI_TJ": dict(min_gc=0.55, min_oe=0.65, min_len=500),
}


def qualifying_window_mask(arr, window, min_gc, min_oe, min_acgt_frac):
    """Boolean mask over window START positions: does the window of `window`
    bases starting here meet the GC / Obs-Exp / non-N criteria?"""
    is_c = (arr == ord("C"))
    is_g = (arr == ord("G"))
    is_acgt = is_c | is_g | (arr == ord("A")) | (arr == ord("T"))
    is_cpg = np.zeros(arr.size, dtype=bool)
    is_cpg[:-1] = is_c[:-1] & is_g[1:]

    def cumsum32(x):
        out = np.zeros(x.size + 1, dtype=np.int32)
        np.cumsum(x, dtype=np.int32, out=out[1:])
        return out

    cum_c, cum_g, cum_cpg, cum_acgt = (cumsum32(v) for v in (is_c, is_g, is_cpg, is_acgt))
    n_starts = arr.size - window + 1
    if n_starts <= 0:
        return np.zeros(0, dtype=bool), None, None, None

    idx = np.arange(n_starts)
    n_c = (cum_c[idx + window] - cum_c[idx]).astype(np.float64)
    n_g = (cum_g[idx + window] - cum_g[idx]).astype(np.float64)
    # dinucleotides starting inside the window end one base before its end
    n_cpg = (cum_cpg[idx + window - 1] - cum_cpg[idx]).astype(np.float64)
    n_acgt = (cum_acgt[idx + window] - cum_acgt[idx]).astype(np.float64)

    gc_frac = (n_c + n_g) / window
    with np.errstate(divide="ignore", invalid="ignore"):
        obs_exp = np.where((n_c > 0) & (n_g > 0), (n_cpg * window) / (n_c * n_g), 0.0)

    mask = (
        (gc_frac >= min_gc)
        & (obs_exp >= min_oe)
        & ((n_acgt / window) >= min_acgt_frac)
    )
    return mask, gc_frac, obs_exp, n_cpg


def merge_runs(mask, window, max_gap=100):
    """Contiguous qualifying window starts -> merged (start, end) regions,
    then gap-merge anything within max_gap of the next region.

    Found 2026-09-11: without gap-merging, a single real CpG-rich region
    fragments into dozens of overlapping near-duplicate "islands" a few bp
    apart - the qualifying mask flickers by 1-2bp right at the threshold
    boundary as the window slides (Obs/Exp oscillating around 0.60), which
    ends one contiguous run and starts another almost immediately. This
    produced 1.2M spurious "islands" genome-wide (vs. the ~20-30K expected
    for a mammalian genome) before this fix. Real implementations (Takai &
    Jones 2002 itself) merge "adjacent or overlapping" qualifying windows -
    a small gap tolerance is standard, not a loosening of the criteria
    themselves (each merged island still has to clear the GC/ObsExp/length
    thresholds on RE-evaluation as a whole region, done by the caller)."""
    if not mask.any():
        return []
    padded = np.concatenate(([False], mask, [False]))
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    raw_regions = []
    for i in range(0, edges.size, 2):
        run_start, run_end = edges[i], edges[i + 1]  # run of qualifying starts
        raw_regions.append([int(run_start), int(run_end - 1 + window)])

    merged = [raw_regions[0]]
    for s, e in raw_regions[1:]:
        if s - merged[-1][1] <= max_gap:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


def region_stats(arr, start, end):
    sub = arr[start:end]
    n_c = int((sub == ord("C")).sum())
    n_g = int((sub == ord("G")).sum())
    is_c = sub == ord("C")
    is_g = sub == ord("G")
    n_cpg = int((is_c[:-1] & is_g[1:]).sum())
    length = end - start
    gc_frac = (n_c + n_g) / length if length else 0.0
    obs_exp = (n_cpg * length) / (n_c * n_g) if n_c > 0 and n_g > 0 else 0.0
    return length, gc_frac, obs_exp, n_cpg


def independent_calls(chrom, per_type_regions):
    """Overlap alone never implies that an interval satisfies both definitions."""
    return [(chrom, s, e, name, length, gc, oe, n)
            for name, regions in per_type_regions.items()
            for s, e, length, gc, oe, n in regions]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--min-acgt-frac", type=float, default=0.9)
    ap.add_argument("--min-contig-size", type=int, default=10_000)
    args = ap.parse_args()

    fa = pysam.FastaFile(args.fasta)
    rows = []
    for chrom, size in zip(fa.references, fa.lengths):
        if size < args.min_contig_size:
            continue
        seq = fa.fetch(chrom).upper()
        arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        del seq

        per_type_regions = {}
        for type_name, p in PARAMS.items():
            mask, _, _, _ = qualifying_window_mask(
                arr, args.window, p["min_gc"], p["min_oe"], args.min_acgt_frac
            )
            regions = merge_runs(mask, args.window)
            # re-evaluate each merged region as a whole, and enforce length
            kept = []
            for s, e in regions:
                length, gc_frac, obs_exp, n_cpg = region_stats(arr, s, e)
                if length >= p["min_len"] and gc_frac >= p["min_gc"] and obs_exp >= p["min_oe"]:
                    kept.append((s, e, length, gc_frac, obs_exp, n_cpg))
            per_type_regions[type_name] = kept
            del mask, regions

        # Preserve each definition's own coordinates and statistics.
        rows.extend(independent_calls(chrom, per_type_regions))

        print(f"{chrom}: {len(per_type_regions['CGI_GGF'])} GGF, "
              f"{len(per_type_regions['CGI_TJ'])} TJ")
        del arr

    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    with open(args.out_bed, "w", encoding="utf-8") as f:
        f.write("#chrom\tstart\tend\ttype\tlength\tgc_frac\tobs_exp\tn_cpg\n")
        for chrom, s, e, type_name, length, gc_frac, obs_exp, n_cpg in rows:
            f.write(f"{chrom}\t{s}\t{e}\t{type_name}\t{length}\t"
                    f"{gc_frac:.4f}\t{obs_exp:.4f}\t{n_cpg}\n")

    n_both = sum(1 for r in rows if r[3] == "CGI_GGF+CGI_TJ")
    n_ggf = sum(1 for r in rows if r[3] == "CGI_GGF")
    n_tj = sum(1 for r in rows if r[3] == "CGI_TJ")
    print(f"\nWrote {len(rows)} CpG islands to {args.out_bed}")
    print(f"  {n_both} called by both parameterizations, {n_ggf} GGF calls, {n_tj} TJ calls (independent)")


if __name__ == "__main__":
    main()
