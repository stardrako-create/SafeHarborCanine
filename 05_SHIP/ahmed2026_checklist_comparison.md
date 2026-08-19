# Checklist comparison against Ahmed et al. 2026 (Cells)

Ahmed, A., Di Molfetta, D., Iaconisi, G.N., et al. (2026). "Human Genome
Safe Harbor Sites: A Comprehensive Review of Criteria, Discovery, Features,
and Applications." *Cells*, 15(1), 81.
DOI: [10.3390/cells15010081](https://doi.org/10.3390/cells15010081)

This review's Figure 1 (Box 1/Box 2 content is embedded as an image, not
extractable text — rendered to read it) lays out 8 core SHS selection
criteria, synthesized from Sadelain et al. 2011 and later refinements
(Pellenz et al.'s "eight SHS criteria" for SHS231, Aznauryan et al.'s Rogi1/
Rogi2 filtering), plus a "Box 1" of additional proposed criteria and a
"Box 2" of acknowledged field-wide challenges.

**We do not currently satisfy all 8 checkboxes.** Honest comparison below.

## The 8 core criteria

| # | Criterion | Our status | Detail |
|---|---|---|---|
| 1 | Distance >=50kb from cancer-unrelated genes | **Yes, reframed** (V3, 2026-08-19) | A literal reading is structurally unsatisfiable for a 50-75kb window that touches genes at both edges by construction (no position inside can be >=50kb from *both* flanking genes at once). Implemented instead as `veto_gene_dense_neighborhood`: is there a THIRD gene within 50kb of either window edge, beyond the two that define it — the real safety concern. Severe: 380/461 candidates excluded (only 45 pass through V3) |
| 2 | Distance >=300kb from cancer-related genes | **Yes** (V2, 2026-08-19) | `veto_risk_gene_radius` — genome-wide 300kb radius search against all 41,632 genes, not just the two flanking genes (199 candidates newly excluded; changed the #1-ranked candidate — see `VERSIONS.md`) |
| 3 | Distance >=300kb from miRNA | **Yes** (2026-08-19) | `scripts/extract_gff3_features.py` pulls all 491 annotated canine miRNA loci from the GFF3; `veto_mirna_nearby` excludes any candidate within 300kb (29 candidates newly excluded) |
| 4 | Outside transcriptional unit | **Yes** | True by construction — SHIP candidates are always intergenic |
| 5 | Outside ultraconserved regions, telomeres, centromeres | Partial | Telomere/centromere GFF3 feature types are excluded during SHIP candidate generation (Ehsan's `features.json`); ultraconserved elements are not checked at all |
| 6 | Outside lncRNA and small RNA | Partial | `ncRNA_gene`-typed features are excluded during SHIP candidate generation at the gene-annotation level; no dedicated lncRNA/small-RNA database cross-check |
| 7 | Located in transcriptionally active (open) chromatin | **Yes** | This is exactly our `moderate_atac` soft-score component (favors accessible-but-not-peak signal) |
| 8 | Outside a TAD containing cancer-related genes | Partial | We veto candidates overlapping a TAD *boundary* and candidates whose flanking genes are risk genes, but we do not check *every* gene inside the candidate's own TAD for cancer relevance — a broader 3D check than what we do |

**Score: 5/8 fully satisfied, 3/8 partial, 0/8 not implemented.** (started at 2/8, 5/8, 1/8; V1 closed criterion 3, V2 closed criterion 2, V3 closed criterion 1 — see `VERSIONS.md`)

## Box 1 — additional proposed criteria

| Criterion | Our status |
|---|---|
| CRISPR/Cas9 editing efficiency and off-target specificity | Not started — explicitly listed as a pending next step (gRNA design phase) in `05_SHIP/README.md` |
| No alteration of transcriptome/proteome/metabolome | Out of scope for this bioinformatic discovery pipeline — belongs to experimental validation (per project scope: "discovery, filtragem, ranking... validação funcional fica a cargo da equipa experimental") |
| No negative impact on stem cell pluripotency/differentiation | Same — experimental validation, not applicable to canine PBMC data anyway |
| "Universal" expression across cell types/tissues | Not tested — our ATAC/RRBS evidence is PBMC-specific; Shrestha et al. 2022 (GEG-SH, already cited in this project) found *zero* shared safe harbors between blood and brain in their own tissue-specific analysis, so "universal" is a genuinely hard bar this pipeline does not claim to clear |

## What this means

Two gaps are worth closing next, roughly in order of how cheap/impactful they are:
1. **miRNA proximity (criterion 3)** — completely missing, and miRNA annotations (`miRNA` feature type) are already present in the RefSeq GFF3 SHIP itself parsed from (`GEF-SH` and this repo's SHIP run excluded `gene`/`ncRNA_gene`/etc. but a dedicated genome-wide miRNA BED + distance check was never built).
2. **300kb radius search (criteria 2 and, in spirit, 8)** — upgrading from "are the two flanking genes risky" to "is any risk gene, or any gene inside the same TAD, within 300kb" is a natural extension of `scripts/score_ship_candidates.py`'s existing risk-gene veto logic.

Ultraconserved elements (part of criterion 5) would need a dedicated
conservation track (e.g. phastCons/phyloP for a multi-species alignment
including dog), which does not currently exist in this project's data.
