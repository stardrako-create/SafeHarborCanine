#!/usr/bin/env python3
"""
Compare two score_ship_candidates_v2.py output tables (e.g. v120b -> v120c)
candidate-by-candidate, keyed on (chrom, start, end) - the fixed 461-
candidate universe doesn't change between versions, only how each one
scores. Written 2026-09-10 so this comparison (already done ad hoc three
times this session: V11->V12->v120b) is a real, reusable script instead of
another scratch one-off.

Reports: candidates that flip hard_veto either direction, which specific
veto(s) and score component(s) caused each flip, and rank/score deltas for
candidates that survive in both versions.
"""
import argparse
import csv


VETO_COLUMNS = [
    "veto_tad_boundary", "veto_atac_peak", "veto_risk_gene", "veto_low_mappability",
    "veto_mirna_nearby", "veto_risk_gene_radius", "veto_gene_dense_neighborhood",
    "veto_lncrna_smallrna", "veto_tad_risk_gene", "veto_high_repeat_content",
    "veto_ultraconserved_element", "veto_external_regulatory_element",
    "veto_low_atac_accessibility",
]
SCORE_COLUMNS = [
    "score_stability_atac", "score_low_peak_frequency", "score_low_methylation",
    "score_stability_rrbs", "score_tad_distance", "final_score", "final_score_percentile",
]


def load(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return {(r["chrom"], r["start"], r["end"]): r for r in rows}


def ranked_survivors(by_key):
    survivors = [(key, row) for key, row in by_key.items() if row["hard_veto"] == "False"]
    survivors.sort(key=lambda kr: -float(kr[1]["final_score"]))
    return {key: i + 1 for i, (key, _) in enumerate(survivors)}


def fmt_locus(key, row):
    genes = f"{row.get('left_gene', '?')}/{row.get('right_gene', '?')}"
    return f"{key[0]}:{key[1]}-{key[2]} ({genes})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--out", default=None, help="optional markdown report path")
    args = ap.parse_args()

    old = load(args.old)
    new = load(args.new)

    old_keys, new_keys = set(old), set(new)
    if old_keys != new_keys:
        only_old = old_keys - new_keys
        only_new = new_keys - old_keys
        print(f"WARNING: candidate universes differ - {len(only_old)} keys only in --old, "
              f"{len(only_new)} only in --new. Expected the same fixed candidate set; "
              f"double-check --old/--new point at the same pipeline run's inputs.")

    shared_keys = old_keys & new_keys
    old_rank = ranked_survivors(old)
    new_rank = ranked_survivors(new)

    newly_passing = [k for k in shared_keys if k not in old_rank and k in new_rank]
    newly_failing = [k for k in shared_keys if k in old_rank and k not in new_rank]
    still_passing = [k for k in shared_keys if k in old_rank and k in new_rank]

    lines = []
    lines.append(f"# Diff: {args.old} -> {args.new}\n")
    lines.append(f"- Old survivors: {len(old_rank)}")
    lines.append(f"- New survivors: {len(new_rank)}")
    lines.append(f"- Newly passing (were vetoed, now survive): {len(newly_passing)}")
    lines.append(f"- Newly failing (used to survive, now vetoed): {len(newly_failing)}")
    lines.append(f"- Survive in both: {len(still_passing)}\n")

    if newly_passing:
        lines.append("## Newly passing")
        for k in sorted(newly_passing, key=lambda k: new_rank[k]):
            old_row, new_row = old[k], new[k]
            failed_vetoes = [v for v in VETO_COLUMNS if old_row.get(v) == "True"]
            lines.append(f"- {fmt_locus(k, new_row)} - new rank #{new_rank[k]}, "
                          f"score={new_row['final_score']}. Previously failed: "
                          f"{', '.join(failed_vetoes) if failed_vetoes else '(none? check hard_veto logic)'}")
        lines.append("")

    if newly_failing:
        lines.append("## Newly failing")
        for k in sorted(newly_failing, key=lambda k: old_rank[k]):
            old_row, new_row = old[k], new[k]
            now_failed = [v for v in VETO_COLUMNS if new_row.get(v) == "True"]
            lines.append(f"- {fmt_locus(k, old_row)} - was rank #{old_rank[k]}, "
                          f"score={old_row['final_score']}. Now fails: "
                          f"{', '.join(now_failed) if now_failed else '(none? check hard_veto logic)'}")
        lines.append("")

    if still_passing:
        lines.append("## Rank/score changes among candidates surviving in both")
        deltas = []
        for k in still_passing:
            old_row, new_row = old[k], new[k]
            rank_delta = old_rank[k] - new_rank[k]  # positive = moved up
            score_delta = float(new_row["final_score"]) - float(old_row["final_score"])
            changed_scores = [
                c for c in SCORE_COLUMNS
                if old_row.get(c) != new_row.get(c)
            ]
            deltas.append((k, old_row, new_row, rank_delta, score_delta, changed_scores))
        deltas.sort(key=lambda d: -abs(d[3]))
        for k, old_row, new_row, rank_delta, score_delta, changed_scores in deltas:
            lines.append(
                f"- {fmt_locus(k, new_row)} - rank #{old_rank[k]} -> #{new_rank[k]} "
                f"({'+' if rank_delta > 0 else ''}{rank_delta}), "
                f"score {old_row['final_score']} -> {new_row['final_score']} "
                f"({'+' if score_delta >= 0 else ''}{score_delta:.4f}), "
                f"changed: {', '.join(changed_scores) if changed_scores else 'none'}"
            )
        lines.append("")

    report = "\n".join(lines)
    print(report)
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"\nWrote report to {args.out}")


if __name__ == "__main__":
    main()
