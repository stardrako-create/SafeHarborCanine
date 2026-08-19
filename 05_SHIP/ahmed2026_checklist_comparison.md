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
| 5 | Outside ultraconserved regions, telomeres, centromeres | Partial | Telomere/centromere GFF3 feature types are excluded during SHIP candidate generation (Ehsan's `features.json`). Ultraconserved elements: investigated 2026-08-19 — no dog-referenced conservation score track (phyloP/phastCons) is readily available. Zoonomia's 241-mammal alignment only publishes projected scores for human (and a few other references); the dog assembly hub on UCSC's CGL server has no conservation bigwig, only the raw 806GB HAL alignment, which would need specialized tooling (e.g. halPhyloP) to project onto dog coordinates — disproportionate effort for this one criterion. Left open rather than forcing a low-quality proxy. |
| 6 | Outside lncRNA and small RNA | **Yes** (V4, 2026-08-19) | `veto_lncrna_smallrna` — direct overlap check against 26,899 lnc_RNA/tRNA/snoRNA/snRNA/guide_RNA/rRNA/SRP_RNA/RNase_P_RNA features. 0 candidates newly excluded — the survivors of V3's gene-density filter were already clean of these |
| 7 | Located in transcriptionally active (open) chromatin | **Yes** | This is exactly our `moderate_atac` soft-score component (favors accessible-but-not-peak signal) |
| 8 | Outside a TAD containing cancer-related genes | **Yes** (V5, 2026-08-19) | `veto_tad_risk_gene` — TAD intervals derived from consecutive boundary calls (`scripts/build_tad_intervals.py`, 6,853 TADs); veto if ANY risk gene falls anywhere inside the candidate's own TAD, not just the two flanking genes or a fixed radius. Only 2 candidates newly excluded — the survivors of V1-V4 already had clean TADs |

**Final score: 7/8 fully satisfied, 1/8 partial (telomere/centromere half of criterion 5 is done; the ultraconserved-elements half stays open — no usable dog-referenced conservation track exists publicly, see row above).** Started at 2/8, 5/8, 1/8; V1 closed criterion 3, V2 closed criterion 2, V3 closed criterion 1, V4 closed criterion 6, V5 closed criterion 8 — see `VERSIONS.md`. This is as far as this pipeline can honestly go against the Ahmed et al. 2026 checklist with data that actually exists for the dog genome today.

## Box 1 — additional proposed criteria

| Criterion | Our status |
|---|---|
| CRISPR/Cas9 editing efficiency and off-target specificity | Not started — explicitly listed as a pending next step (gRNA design phase) in `05_SHIP/README.md` |
| No alteration of transcriptome/proteome/metabolome | Out of scope for this bioinformatic discovery pipeline — belongs to experimental validation (per project scope: "discovery, filtragem, ranking... validação funcional fica a cargo da equipa experimental") |
| No negative impact on stem cell pluripotency/differentiation | Same — experimental validation, not applicable to canine PBMC data anyway |
| "Universal" expression across cell types/tissues | Not tested — our ATAC/RRBS evidence is PBMC-specific; Shrestha et al. 2022 (GEG-SH, already cited in this project) found *zero* shared safe harbors between blood and brain in their own tissue-specific analysis, so "universal" is a genuinely hard bar this pipeline does not claim to clear |

## What's left

Only the ultraconserved-elements half of criterion 5. It would need a
dog-referenced phastCons/phyloP conservation track, which does not exist
publicly today (investigated 2026-08-19 — see row above) — pick this back
up if/when Zoonomia or another project publishes one, or if projecting the
existing human-referenced Zoonomia phyloP through a whole-genome alignment
to dog becomes worth the effort.

Box 1's CRISPR/gRNA criterion remains the actual next phase of this
project (gRNA design and off-target scoring), once a final shortlist is
chosen from `candidates_scored_v5.tsv`.
