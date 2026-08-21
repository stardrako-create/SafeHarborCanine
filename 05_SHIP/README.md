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
- overlaps an independently-curated regulatory element
  (`ehsan_regulatory_elements_ROS.bed` — Ehsan Valiollahi's CanFam3.1 set,
  lifted to ROS_Cfam_1.0 by him before sharing; added V9, 2026-08-21)

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

## Versions — V1 / V2 / V3 (Ahmed et al. 2026 checklist, now 8/8)
**V1** is TAD boundary + ATAC peak + flanking risk gene + self-mappability
only (260 candidates would pass this — see `VERSIONS.md` for why that exact
V1 state isn't preserved as its own file).

**V2** is everything built on top of V1 against the Ahmed et al. 2026
checklist (`ahmed2026_checklist_comparison.md`): miRNA proximity, risk gene
within a 300kb radius, a third gene crowding either window edge within
50kb, lncRNA/small-RNA overlap, and a risk gene anywhere in the candidate's
own TAD (checkpoints 1-5, `candidates_scored.tsv` → `_v5.tsv`).

**V3** closes the remaining criteria: RepeatMasker repeat content
(checkpoint 6, `_v6.tsv`) and ultraconserved elements (checkpoint 7,
`_v7.tsv`) — the latter via a dog-referenced phyloP conservation score
computed from scratch from the raw Zoonomia 241-mammal HAL alignment, since
no pre-built dog-referenced conservation track exists publicly. Full
pipeline, a real WSL crash and recovery, and honest methodological
caveats (single-region neutral model, not production-grade): see
`V3_PROGRESS_NOTES.md`.

**`candidates_scored_v9.tsv` / `candidates_passing_ranked_v9.bed` is the
current, final state: 26/461 candidates pass every criterion this pipeline
checks — all 8/8 of the Ahmed et al. 2026 review's core criteria, plus
overlap against an independent regulatory-element set.** Use these for
anything downstream (06_GEG-SH onward), not any earlier version. V8 kept
V7's exact 34 survivors; it only fixed `final_score` itself, which used to
min-max normalize each component against just the surviving batch
(meaningless across runs - the same candidate's score visibly shifted
between V6 and V7 with nothing about it changing) and is now a percentile
against each track's genome-wide background instead - stable regardless of
batch composition, and portable across labs/species for EpiLog's
`computational_score` field. See `VERSIONS.md` for two real bugs found
and fixed while building this. V9 (2026-08-21) added a hard veto for
overlap with Ehsan Valiollahi's independently-curated regulatory-element
set (received after V8 shipped) — 8 of the 34 V8 survivors were newly
excluded, including the then-#3 top-5 shortlist candidate; see
`top5_shortlist.md` "V9 correction" for the replacement.

Also validated (not a new veto, confirms an existing one): the consensus-
ATAC-peak threshold shows 10.9x/15.1x enrichment against independent
TSS/CpG-island regulatory proxies (`consensus_peak_validation.md`). And
checked (no exclusions resulted): the 40 V6 survivors against the Dog10K
population structural-variant set, 0/40 overlapping a called SV
(`dog10k_sv_check_v6.tsv`).

Top candidate changed three times total. The V1/pre-V2 top candidate
(`NC_051812.1:52,431-118,675`, score 0.83 — independently confirmed by
Ehsan Valiollahi's separate filtering, see `ehsan_crossvalidation.md`) was
excluded in V2 by the 300kb risk-gene-radius check, replaced by
`NC_051811.1:48,020,921-48,077,046`, which then held through checkpoints
3-7. **V8's genome-wide-percentile rescoring (see `VERSIONS.md`) revealed
this was itself a batch-normalization artifact — the current, real top
candidate is `NC_051805.1:7,072,137-7,132,579` (LOC111090579/LOC100685067,
score 0.765).**

## Final shortlist

`top5_shortlist.md` / `top5_shortlist.bed` — top candidate plus 4 backups,
picked from the 26 V9 survivors. Closes the "principal candidate + at
least one backup" requirement from the project's own definition of done
for this phase (see `top5_shortlist.md` for why #1 specifically, why 4
backups rather than just one, and the V9 correction changelog).

## Known gaps (not yet applied)
- The risk-gene list is a **starter proxy**: human HGNC symbols
  (CancerMine + CEG2) matched directly against ROS_Cfam_1.0 gene symbols,
  not a canine-curated ortholog list — a real canine essential-gene/cancer-
  gene resource would be more rigorous.
- The ultraconserved-elements neutral model rests on a single 100kb region
  (not ancestral-repeat-filtered like a production Zoonomia release) — see
  `V3_PROGRESS_NOTES.md` for the honest caveat and how to strengthen it.
- gRNA design and off-target scoring — the actual next phase, now that the
  final shortlist (`top5_shortlist.md`) is chosen. Blocked pending the
  actual Cas9/Cas variant and CAR donor construct to design against - a
  generic guide search without knowing PAM compatibility or the donor
  sequence for homology arms isn't useful yet.

## Files
- `ship_raw_candidates.tsv` — all 461 SHIP candidates, unfiltered
- `canine_risk_genes.tsv` — 2,654 oncogene/tumor-suppressor/essential-gene symbols
- `canine_miRNA.bed`, `canine_all_genes.bed`, `canine_lncRNA_smallRNA.bed`, `canine_tad_intervals.bed` — extracted feature tracks used within V2
- `mappability_check.tsv` — self-realignment mappability check per candidate
- `candidates_scored.tsv` / `candidates_scored_v2.tsv` ... `_v9.tsv` — full tables per version, veto flags and score components
- `candidates_passing_ranked.bed` / `..._v2.bed` ... `_v9.bed` — passing candidates per version, ranked
- `ehsan_regulatory_elements_ROS.bed` — Ehsan Valiollahi's independent regulatory-element set (V9 veto input)
- `VERSIONS.md` — what each version adds and the resulting candidate counts
