# 05_SHIP scoring versions

Each criterion added from `ahmed2026_checklist_comparison.md` gets its own
version — earlier versions' output files are never overwritten, so the
effect of each individual criterion stays auditable.

| Version | Adds | Files | Candidates passing | Top candidate |
|---|---|---|---|---|
| V1 | Baseline (TAD boundary, ATAC peak, flanking risk gene, self-mappability) + miRNA proximity (criterion 3) | `candidates_scored.tsv`, `candidates_passing_ranked.bed` | 243 / 461 | `NC_051812.1:52,431-118,675` (score 0.83) |
| V2 | + risk gene within 300kb radius, not just flanking genes (criterion 2) | `candidates_scored_v2.tsv`, `candidates_passing_ranked_v2.bed` | 164 / 461 | `NC_051817.1:42,409,351-42,467,972` (score 0.75) — **top candidate changed**: the V1 top candidate is now excluded (NLRP3, a cancer/inflammation-associated gene, sits 214kb away — inside the radius, but was never one of its two immediate flanking genes) |
| V3 | + a third gene within 50kb of either window edge (criterion 1, reframed — see comparison doc) | `candidates_scored_v3.tsv`, `candidates_passing_ranked_v3.bed` | 45 / 461 | `NC_051811.1:48,020,921-48,077,046` (score 0.77) — **top candidate changed again**. Severe drop (380 newly excluded) — confirmed with the user (2026-08-19) to keep the hard veto as-is rather than softening it, before continuing to V4-V6 |
| V4 | + overlap with lncRNA/small RNA genes (criterion 6) | `candidates_scored_v4.tsv`, `candidates_passing_ranked_v4.bed` | 45 / 461 | Unchanged from V3 — 0 candidates newly excluded, the V3 survivors were already clean |

See `ahmed2026_checklist_comparison.md` for what each criterion means and
where it comes from.
