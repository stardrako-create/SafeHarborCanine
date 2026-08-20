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
- **V2 (this document) → GitHub release `v0.2.0`, published**
- V3 (ultraconserved elements) → **`v1.0.0`**, in progress

## V3 status (2026-08-20)

Computing a dog-referenced phyloP conservation track from scratch, since none
exists publicly (Zoonomia only distributes human-referenced scores). Pipeline:
download the raw 241-mammal Cactus HAL alignment (806GB) → `hal2maf`
referenced to dog → RepeatMasker on ancestral repeats → `phyloFit` → `phyloP`,
following the Zoonomia consortium's own scripts. HAL download in progress
(resumable, `curl -C -`).

Two more candidate filters were identified and pilot-tested against the
current 43 survivors while V3 downloads, but **deliberately held until V3 is
applied first** — no point spending effort evaluating candidates V3 may
still exclude:

- **RepeatMasker repeat content** (finer-grained than the existing
  self-mappability check, which only catches gross multi-mapping): piloted
  on the 43 current survivors, mean 36.3% repeat content (in line with the
  dog genome average), 3 candidates above 50% (LINE/L1-dominated). Data in
  `repeat_content_v5candidates.tsv`; not yet wired into `score_ship_candidates.py`.
- **Dog10K population structural variants** (Manta-SV, 1,879 dogs,
  `kiddlabshare.med.umich.edu/dog10K/Manta-SV_2022-03-28`): blocked on an
  assembly mismatch — the VCF is UU_Cfam_GSD_1.0, not ROS_Cfam_1.0, with no
  public chain file between them (confirmed via mismatched chromosome
  lengths). Plan: once V3 has reduced the survivor set further, realign just
  those few candidates against UU_Cfam_GSD_1.0 with `bwa mem` (same technique
  already used for the original chr12 candidate) rather than building a full
  genome-wide liftover.

Planned order: **V3 (phyloP) → RepeatMasker → Dog10K SV**, each applied to
the shrinking survivor set from the previous step.
