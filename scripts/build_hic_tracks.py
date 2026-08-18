#!/usr/bin/env python3
"""
Build genome-wide Hi-C structural context tracks from the 3 merged Mischka
Hi-C libraries (PRJNA587469), for use as a "do not disrupt" structural layer
in genomic safe harbor scoring - not a per-dog population signal like the
ATAC/RRBS Mother Tracks, but a single-individual reference map of 3D genome
architecture on ROS_Cfam_1.0 (this Hi-C data predates and structurally
underlies ROS_Cfam_1.0's own chromosome-level scaffolding via RagTag against
Mischka's UU_Cfam_GSD_1.0 assembly).

Pipeline: pairtools merge (3 libraries -> 1) -> cooler cload (1kb fixed bins,
no restriction fragments - the Dovetail kit's enzyme wasn't disclosed in the
source paper) -> cooler zoomify (4DN multi-resolution ladder, ICE-balanced at
every zoom level) -> cooltools insulation (boundary strength + calls, at
multiple window sizes) -> cooltools expected-cis + dots (chromatin loop
anchors, best-effort - skipped with a warning if depth/signal is
insufficient, since this dataset was built for assembly scaffolding, not
dedicated loop-calling depth).

Each step checks whether its output already exists before running, so a
crashed/interrupted run can simply be re-invoked to resume.
"""
import argparse
import os
import subprocess

import pandas as pd
import yaml


def run(cmd):
    print(f"$ {' '.join(str(c) for c in cmd)}", flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def merge_valid_pairs(per_lib_dir, samples, merged_path, nproc, tmpdir):
    if os.path.isfile(merged_path):
        print(f"[skip] {merged_path} already exists")
        return
    inputs = [os.path.join(per_lib_dir, s, "pairs", f"{s}.valid.pairs.gz") for s in samples]
    for p in inputs:
        if not os.path.isfile(p):
            raise SystemExit(f"Missing valid pairs file: {p}")
    os.makedirs(tmpdir, exist_ok=True)
    run(["pairtools", "merge", "--nproc", nproc, "--tmpdir", tmpdir, "-o", merged_path, *inputs])


def build_cool(merged_path, chrom_sizes_path, base_bin, cool_path, assembly):
    if os.path.isfile(cool_path):
        print(f"[skip] {cool_path} already exists")
        return
    run([
        "cooler", "cload", "pairs",
        "--assembly", assembly,
        "-c1", 2, "-p1", 3, "-c2", 4, "-p2", 5,
        f"{chrom_sizes_path}:{base_bin}",
        merged_path,
        cool_path,
    ])


def zoomify_and_balance(cool_path, mcool_path, nproc, resolutions="4DN"):
    if os.path.isfile(mcool_path):
        print(f"[skip] {mcool_path} already exists")
        return
    run([
        "cooler", "zoomify",
        "-r", resolutions,
        "-n", nproc,
        "--balance", "--balance-args", f"-p {nproc}",
        "-o", mcool_path,
        cool_path,
    ])


def build_insulation_view(chrom_sizes_path, windows, view_path):
    """cooltools insulation computes a diagonal diamond score per bin using
    the largest requested window on each side, so a region shorter than
    ~2x the largest window has no valid bins - cooler.annotate() then
    indexes into an empty dataframe and crashes. ROS_Cfam_1.0 has 376
    contigs, most of them tiny unplaced scaffolds far below that
    threshold, so insulation must be restricted to a --view of contigs
    long enough for the window sizes actually requested."""
    min_len = 2 * max(windows)
    sizes = pd.read_csv(chrom_sizes_path, sep="\t", header=None, names=["chrom", "length"])
    kept = sizes[sizes["length"] >= min_len].copy()
    kept["start"] = 0
    kept["name"] = kept["chrom"]
    kept[["chrom", "start", "length", "name"]].to_csv(view_path, sep="\t", header=False, index=False)
    print(f"[view] kept {len(kept)}/{len(sizes)} contigs >= {min_len}bp for insulation")


def call_insulation(mcool_path, insulation_res, windows, nproc, out_tsv, out_dir, chrom_sizes_path):
    if os.path.isfile(out_tsv):
        print(f"[skip] {out_tsv} already exists")
        return
    uri = f"{mcool_path}::resolutions/{insulation_res}"
    view_path = os.path.join(out_dir, "insulation_view.bed")
    build_insulation_view(chrom_sizes_path, windows, view_path)
    cwd = os.getcwd()
    os.chdir(out_dir)
    try:
        run([
            "cooltools", "insulation",
            "-p", nproc,
            "-o", os.path.basename(out_tsv),
            "--threshold", "Li",
            "--view", view_path,
            "--bigwig",
            uri,
            *[str(w) for w in windows],
        ])
    finally:
        os.chdir(cwd)


def export_tad_boundaries(insulation_tsv, windows, out_bed):
    """cooltools insulation flags one is_boundary_<window> column per window
    size; a bin counts as a boundary call here if ANY window size flags it,
    matching the "avoid this position" intent (union, not intersection - we
    want to be conservative about what we treat as structurally sensitive)."""
    df = pd.read_csv(insulation_tsv, sep="\t")
    boundary_cols = [c for c in df.columns if c.startswith("is_boundary_")]
    if not boundary_cols:
        print("[warn] no is_boundary_* columns found in insulation output, skipping BED export")
        return
    is_boundary = df[boundary_cols].fillna(False).any(axis=1)
    boundaries = df.loc[is_boundary, ["chrom", "start", "end"]]
    boundaries.to_csv(out_bed, sep="\t", header=False, index=False)
    print(f"Wrote {len(boundaries)} TAD boundary calls to {out_bed}")


def call_loops(mcool_path, loop_res, nproc, expected_tsv, loops_bedpe, max_loci_sep):
    """Best-effort: this Hi-C data was generated for assembly scaffolding,
    not dedicated deep loop-calling, so cooltools dots may find few or no
    significant dots at the depth available. Failure here should not block
    the rest of the pipeline."""
    if os.path.isfile(loops_bedpe):
        print(f"[skip] {loops_bedpe} already exists")
        return
    uri = f"{mcool_path}::resolutions/{loop_res}"
    try:
        if not os.path.isfile(expected_tsv):
            run(["cooltools", "expected-cis", "-p", nproc, "-o", expected_tsv, uri])
        run([
            "cooltools", "dots",
            "-p", nproc,
            "--max-loci-separation", max_loci_sep,
            "-o", loops_bedpe,
            uri,
            f"{expected_tsv}::balanced.avg",
        ])
    except subprocess.CalledProcessError as e:
        print(f"[warn] loop calling failed or found no significant dots at this depth "
              f"(non-fatal, this track is best-effort): {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="config_hic.yaml with paths/params")
    ap.add_argument("--out-dir", default=None, help="defaults to paths.final_dir in the config")
    ap.add_argument("--base-bin", type=int, default=1000, help="cload base bin size (bp)")
    ap.add_argument("--insulation-res", type=int, default=25000)
    ap.add_argument("--insulation-windows", type=int, nargs="+", default=[100000, 250000, 500000])
    ap.add_argument("--loop-res", type=int, default=10000)
    ap.add_argument("--max-loci-separation", type=int, default=2000000)
    ap.add_argument("--assembly", default="ROS_Cfam_1.0")
    ap.add_argument("--nproc", type=int, default=10)
    args = ap.parse_args()

    with open(args.config, encoding="utf-8") as f:
        config = yaml.safe_load(f)

    samples = config["samples"]
    work_dir = config["paths"]["work_dir"]
    per_lib_dir = os.path.join(work_dir, "per_lib")
    ref_dir = config["paths"]["ref_dir"]
    chrom_sizes_path = os.path.join(ref_dir, "chrom.sizes")
    out_dir = args.out_dir or config["paths"]["final_dir"]
    os.makedirs(out_dir, exist_ok=True)

    merged_path = os.path.join(work_dir, "mischka_merged.valid.pairs.gz")
    cool_path = os.path.join(work_dir, f"mischka.{args.base_bin}.cool")
    mcool_path = os.path.join(out_dir, "mischka_hic.mcool")
    insulation_tsv = os.path.join(out_dir, "tad_insulation.tsv")
    boundaries_bed = os.path.join(out_dir, "tad_boundaries.bed")
    expected_tsv = os.path.join(work_dir, "mischka_expected_cis.tsv")
    loops_bedpe = os.path.join(out_dir, "chromatin_loops.bed")

    print(f"[1/5] Merging {len(samples)} libraries' valid pairs")
    merge_valid_pairs(per_lib_dir, samples, merged_path, args.nproc,
                       tmpdir=os.path.join(work_dir, "tmp_merge"))

    print(f"[2/5] Building base cooler at {args.base_bin}bp")
    build_cool(merged_path, chrom_sizes_path, args.base_bin, cool_path, args.assembly)

    print("[3/5] Zoomifying + balancing (4DN resolution ladder)")
    zoomify_and_balance(cool_path, mcool_path, args.nproc)

    print(f"[4/5] Calling TAD insulation at {args.insulation_res}bp, "
          f"windows={args.insulation_windows}")
    call_insulation(mcool_path, args.insulation_res, args.insulation_windows,
                     args.nproc, insulation_tsv, out_dir, chrom_sizes_path)
    export_tad_boundaries(insulation_tsv, args.insulation_windows, boundaries_bed)

    print(f"[5/5] Calling chromatin loops at {args.loop_res}bp (best-effort)")
    call_loops(mcool_path, args.loop_res, args.nproc, expected_tsv, loops_bedpe,
               args.max_loci_separation)

    print("Done:")
    print(f"  {mcool_path}")
    print(f"  {insulation_tsv}")
    for w in args.insulation_windows:
        print(f"  {os.path.join(out_dir, f'tad_insulation.{w}.bw')}")
    print(f"  {boundaries_bed}")
    if os.path.isfile(loops_bedpe):
        print(f"  {loops_bedpe}")


if __name__ == "__main__":
    main()
