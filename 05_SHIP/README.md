# 05_SHIP — Safe Harbor Integration Prioritization

Per-locus scoring that combines the ATAC-seq, RRBS, and Hi-C tracks from
`04_tracks_processadas/` into a candidate genomic safe harbor shortlist.
Two-stage design, mirroring the two papers this project is built on:

1. **Candidate generation** ([SHIP](https://github.com/MCLeitao/Ship), Leitão
   et al. 2025): intergenic intervals between convergent-orientation genes on
   ROS_Cfam_1.0 (GFF3 accessions `NC_051805.1`–`NC_051844.1`), 50–75 kb.
   **Ehsan Valiollahi (Vasco Barreto lab) ran SHIP itself for the dog
   genome** — filtering the RefSeq GFF3 to the 40 true chromosomes, adapting
   `features.json`, choosing the convergent orientation and 50–75kb size
   range — producing SHIP's actual tool output (`ship_raw_candidates.tsv`,
   461 candidates, parsed here by `scripts/parse_ship_candidates.py`). SHIP's
   own UCSC/Ensembl/regulatory-build cross-referencing was not used in this
   repository — those databases are built for human/mouse/yeast and have no
   dog equivalent of comparable quality; this repository's own ATAC/RRBS/
   Hi-C tracks replace that role entirely (see `ehsan_crossvalidation.md` for
   how this pipeline's filtering compares against Ehsan's own independent
   final shortlist).
2. **Biological filtering and ranking** ([GEG-SH](https://github.com/dewshr/GEG-SH),
   Shrestha et al. 2022 framework, our own data): `scripts/score_ship_candidates.py`
   applies hard vetoes and a transparent soft score using the tracks built in
   this repository.

## Hard vetoes
A candidate is excluded outright if it:
- overlaps a Hi-C TAD boundary (`04_tracks_processadas/.../HiC/tad_boundaries.bed`)
- overlaps an ATAC consensus peak (`.../ATAC/mother_track/consensus_peaks_ATAC.bed`)
  — an unannotated regulatory element sitting inside a nominally intergenic window
- either flanking gene is a known oncogene, tumor suppressor, or core
  essential gene (`canine_risk_genes.tsv`, built by
  `scripts/build_canine_risk_genes.py` from CancerMine — Lever et al. 2019,
  CC0 — and CEG2 — Hart et al. 2017; symbol-matched as a starter proxy for
  canine orthology, see "Known gaps")
- the candidate's own sequence realigns ambiguously against the genome
  (`mappability_check.tsv`, built by `scripts/check_candidate_mappability.py`
  via self-realignment with `bwa mem` — MAPQ < 30 or a secondary/
  supplementary/multi-mapping hit)

## Soft score
Every surviving candidate gets 5 components, each min-max normalized to
[0,1] across the surviving set and averaged (equal weights by default,
tunable via CLI flags — nothing here is a black box, read the formula in
`scripts/score_ship_candidates.py`):
- **stability_atac** — low ATAC signal variability across the 71 dogs
- **stability_rrbs** — low methylation variability across the 71 dogs
- **low_methylation** — low weighted-mean CpG methylation
- **tad_distance** — far from the nearest TAD boundary
- **moderate_atac** — accessibility near the population median (neither
  closed nor unusually open — SHIP already picked intergenic windows, so a
  strong signal here would suggest an unannotated element, not a promoter)

## Versions — V1 vs. V2 (Ahmed et al. 2026 checklist)
**V1** is TAD boundary + ATAC peak + flanking risk gene + self-mappability
only (260 candidates would pass this — see `VERSIONS.md` for why that exact
V1 state isn't preserved as its own file).

**V2** is everything built on top of V1 against the Ahmed et al. 2026
checklist (`ahmed2026_checklist_comparison.md`): miRNA proximity, risk gene
within a 300kb radius, a third gene crowding either window edge within
50kb, lncRNA/small-RNA overlap, and a risk gene anywhere in the candidate's
own TAD. Built one criterion at a time, each checkpoint kept as its own
file (`candidates_scored.tsv` → `_v2.tsv` → `_v3.tsv` → `_v4.tsv` →
`_v5.tsv`) so every individual criterion's effect stays auditable — see
`VERSIONS.md` for the full checkpoint-by-checkpoint table.

**`candidates_scored_v5.tsv` / `candidates_passing_ranked_v5.bed` is V2's
final, current state: 43/461 candidates pass every criterion this pipeline
checks.** Use these for anything downstream (06_GEG-SH onward), not the
plain-named files. The pipeline now satisfies 7/8 of the review's core
criteria (only the ultraconserved-elements half of criterion 5 remains
open — no usable dog-referenced conservation track exists publicly,
investigated and documented in the comparison doc).

Top candidate changed twice during V2 (see `VERSIONS.md`): the V1/pre-V2
top candidate (`NC_051812.1:52,431-118,675`, score 0.83 — independently
confirmed by Ehsan Valiollahi's separate filtering, see
`ehsan_crossvalidation.md`) is excluded in V2 by the 300kb risk-gene-radius
check. **V2's top candidate is `NC_051811.1:48,020,921-48,077,046`
(score 0.77).**

## Known gaps (not yet applied)
- The risk-gene list is a **starter proxy**: human HGNC symbols
  (CancerMine + CEG2) matched directly against ROS_Cfam_1.0 gene symbols,
  not a canine-curated ortholog list — a real canine essential-gene/cancer-
  gene resource would be more rigorous.
- The mappability check only catches gross self-multi-mapping (whole
  candidate realigns ambiguously) — it does not measure partial repeat
  *content* within a candidate the way a full RepeatMasker run would; only
  2/461 candidates were flagged, so finer-grained repeat content inside an
  otherwise-unique window is not yet screened.
- Ultraconserved elements (criterion 5) — no dog-referenced conservation
  track exists publicly; see `ahmed2026_checklist_comparison.md`.
- Population variant data (Dog10K) for structural instability at candidates
- gRNA design and off-target scoring — the actual next phase, once a final
  shortlist is chosen from V2's final checkpoint (`candidates_scored_v5.tsv`)

## Files
- `ship_raw_candidates.tsv` — all 461 SHIP candidates, unfiltered
- `canine_risk_genes.tsv` — 2,654 oncogene/tumor-suppressor/essential-gene symbols
- `canine_miRNA.bed`, `canine_all_genes.bed`, `canine_lncRNA_smallRNA.bed`, `canine_tad_intervals.bed` — extracted feature tracks used within V2
- `mappability_check.tsv` — self-realignment mappability check per candidate
- `candidates_scored.tsv` / `candidates_scored_v2.tsv` ... `_v5.tsv` — full tables per version, veto flags and score components
- `candidates_passing_ranked.bed` / `..._v2.bed` ... `_v5.bed` — passing candidates per version, ranked
- `VERSIONS.md` — what each version adds and the resulting candidate counts
