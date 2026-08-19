# 05_SHIP scoring versions

Each criterion added from `ahmed2026_checklist_comparison.md` gets its own
version — earlier versions' output files are never overwritten, so the
effect of each individual criterion stays auditable.

| Version | Adds | Files | Candidates passing | Top candidate |
|---|---|---|---|---|
| V1 | Baseline (TAD boundary, ATAC peak, flanking risk gene, self-mappability) + miRNA proximity (criterion 3) | `candidates_scored.tsv`, `candidates_passing_ranked.bed` | 243 / 461 | `NC_051812.1:52,431-118,675` (score 0.83) |
| V2 | + risk gene within 300kb radius, not just flanking genes (criterion 2) | `candidates_scored_v2.tsv`, `candidates_passing_ranked_v2.bed` | 164 / 461 | `NC_051817.1:42,409,351-42,467,972` (score 0.75) — **top candidate changed**: the V1 top candidate is now excluded (NLRP3, a cancer/inflammation-associated gene, sits 214kb away — inside the radius, but was never one of its two immediate flanking genes) |

See `ahmed2026_checklist_comparison.md` for what each criterion means and
where it comes from.
