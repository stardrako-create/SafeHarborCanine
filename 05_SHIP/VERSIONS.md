# 05_SHIP scoring versions

**V1** = the original scoring (TAD boundary, ATAC peak, flanking risk gene,
self-mappability only) — no Ahmed et al. 2026 checklist criteria applied.

**V2** = everything done against the Ahmed et al. 2026 checklist
(`ahmed2026_checklist_comparison.md`), built up one criterion at a time.
The `_v2`/`_v3`/`_v4`/`_v5` filename suffixes below are the intermediate
checkpoints *within* V2, kept as separate files (never overwriting an
earlier checkpoint) so each individual criterion's effect stays auditable
— they are not separate top-level versions themselves. **V2's final state
is `candidates_scored_v5.tsv` / `candidates_passing_ranked_v5.bed`.**

One naming wrinkle, for transparency: the miRNA criterion (checkpoint 1) was
implemented and written in place to `candidates_scored.tsv` /
`candidates_passing_ranked.bed` *before* the "keep every step as its own
file" instruction was given, so the original pre-Ahmed-checklist V1 numbers
(TAD+ATAC+risk-gene+mappability only, 260 passing) were not preserved as a
separate file — `candidates_scored.tsv` already contains the miRNA
criterion. Everything from checkpoint 2 onward (`_v2.tsv` and later) is
untouched once written.

| Checkpoint (within V2) | Adds | Files | Candidates passing | Top candidate |
|---|---|---|---|---|
| *(pre-V2 baseline, not separately preserved — see note above)* | TAD boundary, ATAC peak, flanking risk gene, self-mappability | — | 260 / 461 | `NC_051812.1:52,431-118,675` (score 0.83) |
| 1: miRNA | + miRNA proximity, 300kb (criterion 3) | `candidates_scored.tsv`, `candidates_passing_ranked.bed` | 243 / 461 | `NC_051812.1:52,431-118,675` (score 0.83) — unchanged |
| 2: risk-gene radius | + risk gene within 300kb radius, not just flanking genes (criterion 2) | `candidates_scored_v2.tsv`, `candidates_passing_ranked_v2.bed` | 164 / 461 | `NC_051817.1:42,409,351-42,467,972` (score 0.75) — **top candidate changed**: NLRP3 (cancer/inflammation-associated) sits 214kb away — inside the radius, but was never one of its two immediate flanking genes |
| 3: gene-dense neighborhood | + a third gene within 50kb of either window edge (criterion 1, reframed — see comparison doc) | `candidates_scored_v3.tsv`, `candidates_passing_ranked_v3.bed` | 45 / 461 | `NC_051811.1:48,020,921-48,077,046` (score 0.77) — **top candidate changed again**. Severe drop (380 newly excluded) — confirmed with the user (2026-08-19) to keep the hard veto as-is rather than softening it |
| 4: lncRNA/smallRNA | + overlap with lncRNA/small RNA genes (criterion 6) | `candidates_scored_v4.tsv`, `candidates_passing_ranked_v4.bed` | 45 / 461 | Unchanged — 0 candidates newly excluded, checkpoint 3's survivors were already clean |
| 5: TAD content | + any risk gene inside the candidate's own TAD, not just flanking/radius (criterion 8) | `candidates_scored_v5.tsv`, `candidates_passing_ranked_v5.bed` | **43 / 461 — V2 final** | Unchanged — only 2 newly excluded, top candidate's TAD is clean |

Ultraconserved elements (part of criterion 5) were investigated but not
implemented — no dog-referenced conservation track exists publicly. See
`ahmed2026_checklist_comparison.md`.

**Bottom line: V2's `candidates_scored_v5.tsv` / `candidates_passing_ranked_v5.bed`
is the current, most rigorous scoring — use these, not the plain-named V1 files,
for anything downstream (06_GEG-SH onward).**

## Release tagging

- V1 (pre-checklist baseline) → no dedicated tag, superseded before release tagging started
- **V2 (this document) → GitHub release `v0.2.0`**
- V3 (ultraconserved elements, if a dog-referenced conservation track becomes
  available — see `ahmed2026_checklist_comparison.md`) → **`v1.0.0`**
