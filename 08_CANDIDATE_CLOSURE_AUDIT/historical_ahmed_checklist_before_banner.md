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
| 1 | Distance >=50kb from cancer-unrelated genes | **Yes, corrected 2026-09-11** (was wrong 2026-08-19 → 2026-09-10) | The paper's actual body text (Section 2, page 3 — verified directly against the PDF, not just this table's earlier paraphrase) is specific: "at least 50 kilobases (kb), from the **5' end** of coding genes" — not gene-body distance. V3's implementation (2026-08-19) measured distance to the nearest gene-body edge, excluding the two flanking genes by name — a different criterion, not a conservative version of the literal one. That exclusion's own justification (a 50-75kb window can never be 50kb from both flanking genes) is true for body-distance but not for 5'-end distance: SHIP only returns convergent gene pairs, so the flanking genes' near edges (facing the window) are their 3' ends by construction — their 5' ends face outward and aren't automatically clear. **Fixed 2026-09-11**: `scripts/extract_genes_stranded.py` extracts protein-coding genes with real strand from the GFF3 (20,950 genes); `veto_gene_dense_neighborhood` now measures distance to the nearest 5' end (strand-aware), among ALL genes, no exclusion — see `distance_to_nearest_5prime_end()` in `score_ship_candidates_v2.py`. Effect (holding every other veto/threshold fixed): gene-density failures 380→264/461. The old body-distance metric is kept as `gene_dense_clearance_old_bodydist` for comparison only, not used in the veto. |
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

## Project-specific additions beyond Ahmed's 8 criteria (not to be described as "Ahmed criteria")

Two hard vetoes in `score_ship_candidates_v2.py` are this project's own
additions, not implementations of anything in Ahmed et al. 2026's checklist
— worth being explicit about this in any future write-up so "passes every
Ahmed criterion" and "passes every implemented veto" aren't conflated
(a distinction a 2026-09-11 review correctly pushed on):

- `veto_atac_peak` / `score_low_peak_frequency` — excludes direct overlap
  with a called ATAC peak, and continuously rewards low peak-frequency
  nearby. Ahmed's criterion 7 is about being IN active/open chromatin
  generally, not a rule against sitting near a discrete called peak — this
  is our own, additional safety margin against disrupting a real
  regulatory element, layered on top of criterion 7, not a restatement
  of it.
- `veto_external_regulatory_element` — hard excludes overlap with
  `ehsan_regulatory_elements_ROS.bed` (see below).

**The regulatory-element file itself was audited 2026-09-11 and found
unauditable as received**: 75,600 bare intervals (chrom/start/end only) —
no element type, no evidence/score, no name, no header. Traced its
provenance: Ehsan's own CanFam3.1 set, lifted to ROS_Cfam_1.0 by him,
received 2026-08-21 (`ehsan_crossvalidation.md`, "Update 2026-08-21"). That
file predates a 2026-09-11 email in which Ehsan cited a regulatory element
at a locus this file does NOT flag — consistent with his own curation
having moved on since, not with an error on our side.

**Since the file can't be audited or currently reconciled, an independent,
self-auditable component was built rather than continuing to depend on it
blind**: `scripts/build_cpg_islands.py` computes CpG islands natively on
the ROS_Cfam_1.0 FASTA (no liftover, so no liftover error), both the
classic Gardiner-Garden & Frommer 1987 (GC>=50%, ObsExp>=0.60, >=200bp) and
the stricter Takai & Jones 2002 (GC>=55%, ObsExp>=0.65, >=500bp)
parameterizations, each region typed and carrying its own GC%/ObsExp/CpG-
count evidence — output: `05_SHIP/cpg_islands_ROS_Cfam_1.0.bed`
(113,459 islands: 74,242 GGF-only, 38,691 satisfying both, 526 TJ-only).
CpG islands were part of Ehsan's own stated methodology too ("CanFam3.1
regulatory elements/**CpG islands** liftover"), so this is comparable to
his approach, computed independently, not an unrelated substitute.

*(Two real bugs were found and fixed while building this, both confirmed
via direct inspection of the output before trusting any count: a 1bp
sliding-window mask that flickered right at the threshold boundary
fragmented single real islands into dozens of overlapping near-duplicates
- 1.2M spurious "islands" before a gap-tolerant merge (100bp) fixed it down
to the biologically plausible 113,459; and the "called by both
parameterizations" check compared exact (start,end) tuples instead of
genomic overlap, undercounting agreement between the two independently-
merged region sets by roughly 1000x - 17 vs. the correct 38,691.)*

**Result: 55.0% of Ehsan's 75,600 elements overlap our independently-built
CpG islands** — well above chance (CpG islands cover only ~1-2% of the
genome), a real, independent validation that a majority of his set is
genuinely promoter/CpG-associated. This is NOT proof his whole set is
correct (CpG islands are a promoter-specific proxy — they don't cover
distal enhancers, likely a large share of the remaining 45%), and it does
NOT resolve the specific rank-#8 dispute: our CpG islands show zero overlap
there, same as his own file — an honest inconclusive (that locus isn't a
CpG island in either dataset, consistent with the disputed element being a
non-CpG-island type our check was never going to catch either way), not a
refutation of his claim. Still needs him to specify element type/evidence
for that one.

Not wired into `hard_veto` as a replacement for `veto_external_regulatory_
element` — it's a complementary, independently-auditable check, not a full
substitute for an enhancer-inclusive annotation. Whether/how to combine the
two is an open design question, not yet decided.

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
