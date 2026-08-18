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
461 SHIP candidates → 181 hard-vetoed (153 TAD boundary, 32 ATAC peak
overlap) → **280 ranked candidates** in `candidates_scored.tsv` (full table)
and `candidates_passing_ranked.bed` (passing only, sorted by score).

## Known gaps (not yet applied)
- RepeatMasker / segmental duplication / mappability filtering
- Canine oncogene / tumor suppressor / essential-gene proximity list
- Population variant data (Dog10K) for structural instability at candidates
- gRNA design and off-target scoring — only after the above land

## Files
- `ship_raw_candidates.tsv` — all 461 SHIP candidates, unfiltered
- `candidates_scored.tsv` — full table with veto flags and score components
- `candidates_passing_ranked.bed` — 280 passing candidates, ranked
