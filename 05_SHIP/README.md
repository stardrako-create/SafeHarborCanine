# 05_SHIP — Safe Harbor Integration Prioritization

Per-locus scoring that combines the ATAC-seq, RRBS, and Hi-C tracks from
`04_tracks_processadas/` into a candidate genomic safe harbor shortlist.
Two-stage design, mirroring the two papers this project is built on:

1. **Candidate generation** ([SHIP](https://github.com/MCLeitao/Ship), Leitão
   et al. 2025): intergenic intervals between convergent-orientation genes on
   ROS_Cfam_1.0 (GFF3 accessions `NC_051805.1`–`NC_051844.1`), 50–75 kb, from
   SHIP's actual tool output (`ship_raw_candidates.tsv`, 461 candidates,
   parsed by `scripts/parse_ship_candidates.py`). SHIP's own UCSC/Ensembl/
   regulatory-build cross-referencing was not used — those databases are
   built for human/mouse/yeast and have no dog equivalent of comparable
   quality; our own ATAC/RRBS/Hi-C tracks replace that role entirely.
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

## Current result
461 SHIP candidates → 201 hard-vetoed (153 TAD boundary, 32 ATAC peak
overlap, 47 risk gene, 2 low mappability — categories overlap, a candidate
can trigger more than one) → **260 ranked candidates** in
`candidates_scored.tsv` (full table) and `candidates_passing_ranked.bed`
(passing only, sorted by score). Top candidate (`NC_051812.1:52,431-118,675`,
score 0.83) is unchanged from before the risk-gene/mappability vetoes were
added, and is independently confirmed by Ehsan Valiollahi's separate
filtering approach — see `ehsan_crossvalidation.md`.

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
- Population variant data (Dog10K) for structural instability at candidates
- gRNA design and off-target scoring — only after the above land

## Files
- `ship_raw_candidates.tsv` — all 461 SHIP candidates, unfiltered
- `canine_risk_genes.tsv` — 2,654 oncogene/tumor-suppressor/essential-gene symbols
- `mappability_check.tsv` — self-realignment mappability check per candidate
- `candidates_scored.tsv` — full table with veto flags and score components
- `candidates_passing_ranked.bed` — 260 passing candidates, ranked
