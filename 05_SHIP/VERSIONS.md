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

| 6: repeat content | + RepeatMasker %% repeat content > 50%% (finer than self-mappability, catches partial repeat content within an otherwise-unique window) | `candidates_scored_v6.tsv`, `candidates_passing_ranked_v6.bed` | 40 / 461 | Unchanged — 3 newly excluded, top candidate's window is clean |
| 7: ultraconserved elements | + phyloP 50bp-rolling-mean conservation > 6.5, dog-referenced, computed from the raw Zoonomia HAL (criterion 5's other half — see `V3_PROGRESS_NOTES.md`) | `candidates_scored_v7.tsv`, `candidates_passing_ranked_v7.bed` | **34 / 461 — current final** | Unchanged coordinates (`NC_051811.1:48,020,921-48,077,046`); soft score shifted 0.774→0.704 because min-max normalization is relative to the current survivor set, not a real change in the underlying tracks |

**All 8/8 Ahmed et al. 2026 checklist criteria now satisfied** (checkpoint 7
closes the last one). See `V3_PROGRESS_NOTES.md` for the full ultraconserved-
elements pipeline and its honest limitations (single-region neutral model,
not production-grade Zoonomia methodology).

**Bottom line: `candidates_scored_v8.tsv` / `candidates_passing_ranked_v8.bed`
is the current, most rigorous scoring — use these, not any earlier version,
for anything downstream (06_GEG-SH onward). Same 34 candidates as V7; V8
only fixed how `final_score` is computed (genome-wide percentile, not
batch min-max — see below), which changed the top candidate.**

## Consensus ATAC peak threshold — validated (2026-08-21)

Ehsan raised a fair methodological question: the consensus-peak veto uses a
≥36/71 (majority) reproducibility threshold following ArchR's convention,
but that convention was never checked against independent evidence that it
actually identifies functional regulatory elements in this dataset.
Validated: intersected `consensus_peaks_ATAC.bed` against TSS±2kb (GFF3) and
CpG islands (computed natively on ROS_Cfam_1.0, not cross-species liftover)
— 10.9x and 15.1x enrichment respectively over genome background. Real
support for the threshold; kept as-is. See `consensus_peak_validation.md`.

## Dog10K structural variants — checked, no exclusions (2026-08-21)

`bwa mem` realignment of the (then-)40 V6 survivors against UU_Cfam_GSD_1.0
(approximate/best-effort, no chain file exists to that assembly — see
`V3_PROGRESS_NOTES.md`), checked against the Dog10K population SV set
(Manta-SV, 1,879 dogs). 40/40 realigned confidently (MAPQ≥30), **0/40
overlapped any called structural variant** — no exclusions, not wired into
`score_ship_candidates.py` as a formal veto input since it excluded nobody.
See `dog10k_sv_check_v6.tsv`.

## Soft-score normalization rewrite — genome-wide percentile, not batch min-max (2026-08-21)

`final_score` used to min-max normalize each soft-score component against
just the surviving candidate batch. That made the number meaningless
across runs: the same candidate's score visibly shifted (0.774 -> 0.704)
between V6 and V7 purely because the survivor pool shrank, with nothing
about the candidate itself changing. It also wasn't portable - a
prerequisite for asking other labs/species to submit comparable
`computational_score` values into EpiLog (epilogbio.netlify.app).

Rewrote to express each component as a percentile rank against that
track's **genome-wide background distribution** instead - stable
regardless of batch composition, and every submitter already has this
denominator (the track itself, pre-filtering). `candidates_scored_v8.tsv`
/ `candidates_passing_ranked_v8.bed` - same 34 survivors as V7 (no vetoes
changed), only the ranking/score changed.

**The top candidate changed as a direct result**: `NC_051811.1:48,020,921-
48,077,046` had been #1 since V2 checkpoint 3, but that was an artifact of
comparing it only to its 33-42 podium-mates. Against the true genome-wide
background, **`NC_051805.1:7,072,137-7,132,579` (LOC111090579/
LOC100685067) is the real top candidate**, score 0.765.

Two real bugs found and fixed while building this, both worth knowing
about if this code is touched again:
1. **Memory**: `pyBigWig`'s `.values()` returns one entry per *base pair*,
   not per stored bin - building a genome-wide background this way meant
   ~2.4 billion points per track instead of the ~48-96 million the tracks
   actually contain, and drove one run to 28GB RAM before it was killed
   (the same failure mode that crashed WSL earlier that night, caught in
   time this time). Fixed by using `.intervals()` instead, which returns
   the track's own native bins.
2. **RRBS background dominated by no-coverage zeros**: RRBS is sparse
   (MspI-site-concentrated) - 88.6% of the genome has zero dogs with
   confident coverage there. An unfiltered background was ~93%
   exactly-zero, so every candidate (which has real coverage, or is
   already excluded via `no_rrbs_coverage`) looked artificially extreme
   against a background that was mostly measuring "not sequenced here,"
   not "genuinely low here." Caught by a sanity check: `score_stability_rrbs`
   was suspiciously near-identical (~0.063) across every top candidate
   before the fix - not real differentiation. Fixed by masking the RRBS
   mean/variability backgrounds to coverage>0 bins only (using the same
   `RRBS_cpg_coverage_frequency.bw` already used per-candidate).

## Release tagging

- V1 (pre-checklist baseline) → no dedicated tag, superseded before release tagging started
- **V2 (checkpoints 1-5) → GitHub release `v0.2.0`, published**
- **V3 (checkpoints 6-7, this document) → `v1.0.0`** — all 8/8 criteria now
  closed; tagging is the user's own action, not done automatically.
