#!/usr/bin/env python3
"""
Scoring-level robustness check for score_ship_candidates_v2.py, written
2026-09-10, REDESIGNED 2026-09-11 after a reviewer correctly caught that the
first version's single pooled "survive %/top-1 %" statistic was misleading:
hard_veto (including veto_low_atac_accessibility) is computed entirely
before weights are ever touched in score_ship_candidates_v2.py's main() -
scoring weights change the ranking among candidates that already passed
every veto, they cannot change WHICH candidates pass. At threshold=0.55,
where only one candidate (ANO2/NTF3) clears the accessibility veto, EVERY
weight-perturbation config is therefore structurally guaranteed to report
that same candidate as the sole survivor and rank #1 - not because it's
robust, but because there was nothing else in the field to compare it
against. Mixing those configs into one "robustness" percentage alongside
genuine threshold variation manufactured an inflated, partly circular
number (verified directly against the code: hard_veto is set in the first
per-candidate loop, final_score/weights only in the second).

Fixed by separating two genuinely different questions instead of pooling
them into one statistic:

1. THRESHOLD SENSITIVITY (weight-independent by construction, since vetoes
   don't depend on weights) - how does the survivor SET change across the
   threshold grid, holding weights at baseline. This is the same question
   the earlier version's "threshold-only" configs asked; now reported on
   its own, not blended with weight-perturbation results.

2. RANKING STABILITY UNDER WEIGHT PERTURBATION (threshold-independent by
   design) - held at ONE fixed, chosen threshold with a real field of
   multiple survivors (not p55, where there's only one candidate and
   nothing to reorder), run weight perturbations and report how much the
   ranking among that FIXED survivor set reshuffles - mean rank, rank
   range, and how often each candidate is #1, among candidates that
   actually had competition.

Every weight and the accessibility threshold are already exposed as CLI
flags on score_ship_candidates_v2.py (see run_score_v120.sh) - this script
drives that same, unmodified entrypoint as a subprocess, so both checks
exercise the literal code path the "official" run uses, not a
reimplementation of the scoring logic.

Usage:
    python3 sensitivity_analysis.py \\
      --score-script score_ship_candidates_v2.py \\
      --fixed-args "--candidates ... (everything run_score_v120.sh passes
        except --min-atac-accessibility-percentile, --w-*, --out-scored,
        --out-passing-bed)" \\
      --out-dir ../06_v1.20.0/sensitivity \\
      --ranking-threshold 0.25
"""
import argparse
import csv
import os
import random
import shlex
import subprocess
import sys

BASE_WEIGHTS = {
    "stability_atac": 1.0,
    "stability_rrbs": 1.0,
    "low_methylation": 1.0,
    "tad_distance": 1.0,
    "low_peak_frequency": 1.0,
}
BASE_THRESHOLD = 0.55
# Range actually under discussion 2026-09-11 (p10-p55) - the accessibility
# recalibration made p55 itself the strict/edge case, not the middle of the
# interesting range.
THRESHOLD_GRID = [0.10, 0.15, 0.20, 0.25, 0.35, 0.45, 0.55]


def run_one_config(score_script, fixed_args, weights, threshold, out_dir, config_name):
    out_scored = os.path.join(out_dir, f"scored_{config_name}.tsv")
    out_bed = os.path.join(out_dir, f"passing_{config_name}.bed")
    cmd = [sys.executable, score_script] + fixed_args + [
        "--min-atac-accessibility-percentile", str(threshold),
        "--w-stability-atac", str(weights["stability_atac"]),
        "--w-stability-rrbs", str(weights["stability_rrbs"]),
        "--w-low-methylation", str(weights["low_methylation"]),
        "--w-tad-distance", str(weights["tad_distance"]),
        "--w-low-peak-frequency", str(weights["low_peak_frequency"]),
        "--out-scored", out_scored,
        "--out-passing-bed", out_bed,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"config {config_name} failed:\n{result.stderr[-2000:]}")
    survivors_path = out_scored
    with open(survivors_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    os.remove(out_scored)
    if os.path.exists(out_bed):
        os.remove(out_bed)
    survivors = [r for r in rows if r["hard_veto"] == "False"]
    survivors.sort(key=lambda r: -float(r["final_score"]))
    return survivors


def threshold_sensitivity(score_script, fixed_args, out_dir):
    """Weight-independent: survivor set/count at each threshold, baseline
    weights throughout. Answers "how does the candidate SET change with the
    threshold" - a question weight perturbation cannot inform, since vetoes
    don't read the weights at all."""
    print("=== 1. Threshold sensitivity (baseline weights, survivor SET per threshold) ===")
    rows_out = []
    for t in THRESHOLD_GRID:
        survivors = run_one_config(score_script, fixed_args, BASE_WEIGHTS, t, out_dir, f"thresh_{t}")
        names = ", ".join(f"{r['chrom']}:{r['start']}-{r['end']} ({r['left_gene']}/{r['right_gene']})"
                           for r in survivors[:5])
        more = f" (+{len(survivors)-5} more)" if len(survivors) > 5 else ""
        print(f"  p{int(t*100)}: {len(survivors)} survive - top by score: {names}{more}")
        rows_out.append((t, survivors))
    return rows_out


def ranking_stability(score_script, fixed_args, out_dir, threshold, perturbation, n_random, seed):
    """Threshold-independent: fix the threshold (and therefore the survivor
    SET) and vary only the weights, to see whether the ranking WITHIN that
    fixed set is stable. Only meaningful with >1 survivor - the caller picks
    a threshold with a real field."""
    print(f"\n=== 2. Ranking stability under weight perturbation, held at p{int(threshold*100)} ===")
    configs = [("baseline", dict(BASE_WEIGHTS))]
    for name in BASE_WEIGHTS:
        for direction, mult in [("minus", 1 - perturbation), ("plus", 1 + perturbation)]:
            weights = dict(BASE_WEIGHTS)
            weights[name] = round(weights[name] * mult, 4)
            configs.append((f"{name}_{direction}{int(perturbation*100)}pct", weights))
    rng = random.Random(seed)
    for i in range(n_random):
        weights = {name: round(v * rng.uniform(1 - perturbation, 1 + perturbation), 4)
                   for name, v in BASE_WEIGHTS.items()}
        configs.append((f"random_{i}", weights))

    ranks = {}
    candidate_info = {}
    survivor_set_sizes = set()
    for i, (config_name, weights) in enumerate(configs, 1):
        print(f"  [{i}/{len(configs)}] {config_name}")
        survivors = run_one_config(score_script, fixed_args, weights, threshold, out_dir,
                                    f"rank_{config_name}")
        survivor_set_sizes.add(len(survivors))
        for rank, r in enumerate(survivors, 1):
            key = (r["chrom"], r["start"], r["end"])
            ranks.setdefault(key, []).append(rank)
            candidate_info.setdefault(key, f"{r['left_gene']}/{r['right_gene']}")

    if len(survivor_set_sizes) > 1:
        print(f"  [warn] survivor set size varied across configs ({survivor_set_sizes}) - "
              f"weights should never change hard_veto; investigate before trusting this report")

    return ranks, candidate_info, len(configs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-script", default="score_ship_candidates_v2.py")
    ap.add_argument("--fixed-args", required=True,
                     help="everything score_ship_candidates_v2.py needs EXCEPT "
                          "--min-atac-accessibility-percentile, --w-*, --out-scored, "
                          "--out-passing-bed - as one shell-quoted string")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ranking-threshold", type=float, default=0.25,
                     help="fixed threshold for the ranking-stability check - pick one with "
                          "a real field of survivors, not p55 (only 1 survivor there)")
    ap.add_argument("--perturbation", type=float, default=0.30)
    ap.add_argument("--n-random", type=int, default=15)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--report", default=None, help="defaults to <out-dir>/sensitivity_report.md")
    args = ap.parse_args()

    fixed_args = shlex.split(args.fixed_args)
    os.makedirs(args.out_dir, exist_ok=True)
    report_path = args.report or os.path.join(args.out_dir, "sensitivity_report.md")

    threshold_rows = threshold_sensitivity(args.score_script, fixed_args, args.out_dir)
    ranks, candidate_info, n_rank_configs = ranking_stability(
        args.score_script, fixed_args, args.out_dir, args.ranking_threshold,
        args.perturbation, args.n_random, args.seed,
    )

    lines = ["# Sensitivity report (redesigned 2026-09-11 - see script docstring for why)\n"]

    lines.append("## 1. Threshold sensitivity (survivor SET per threshold, baseline weights)\n")
    lines.append("| threshold | # survivors |")
    lines.append("|---|---|")
    for t, survivors in threshold_rows:
        lines.append(f"| p{int(t*100)} | {len(survivors)} |")
    lines.append("")

    lines.append(f"## 2. Ranking stability under weight perturbation, held at "
                 f"p{int(args.ranking_threshold*100)} ({n_rank_configs} weight configs)\n")
    lines.append(
        "Threshold is fixed here, so the survivor SET is fixed too (weights "
        "cannot change hard_veto) - this isolates whether the SCORE ORDERING "
        "among that fixed field is stable under weight uncertainty, which "
        "p55's single-survivor field could not test.\n"
    )
    lines.append("| locus | genes | mean rank | rank range | #1 in N/total configs |")
    lines.append("|---|---|---|---|---|")
    for key in sorted(ranks, key=lambda k: sum(ranks[k]) / len(ranks[k])):
        r = ranks[key]
        mean_rank = sum(r) / len(r)
        n1 = sum(1 for x in r if x == 1)
        lines.append(f"| {key[0]}:{key[1]}-{key[2]} | {candidate_info[key]} | "
                      f"{mean_rank:.1f} | {min(r)}-{max(r)} | {n1}/{n_rank_configs} |")

    report = "\n".join(lines)
    print("\n" + report)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"\nWrote {report_path}")


if __name__ == "__main__":
    main()
