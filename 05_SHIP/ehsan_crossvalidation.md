# Cross-validation against Ehsan Valiollahi's independent SHIP filtering

Ehsan Valiollahi (Vasco Barreto lab) independently ran the SHIP tool against
the same ROS_Cfam_1.0 GFF3 (40 true chromosomes, convergent orientation,
50–75 kb intergenic windows → the same 461 candidates this repository's
`ship_raw_candidates.tsv` is parsed from), then applied a different final
filter: mapping known human regulatory regions (promoters, conserved miRNA
targets, TFBS) onto the dog genome via liftover, and excluding any
candidate overlapping one. That narrowed 461 → **27 final candidates**
(`all_safe_harbors_complete_records.txt`, shared 2026-05-04).

This repository's `score_ship_candidates.py` instead uses our own canine
population data (71-dog ATAC/RRBS, single-individual Hi-C), plus a canine
risk-gene proximity check and a self-mappability check, rather than
human-liftover homology, narrowing the same 461 → 260 ranked candidates.

## Result (V1-era, superseded — see V2 update below)

- All 27 of Ehsan's candidates are confirmed present in the shared 461-candidate set.
- **Our #1-ranked candidate (`NC_051812.1:52431-118675`, score 0.83) is also in Ehsan's final 27** — two independent filtering strategies converge on the same top locus, and it stays #1 even after adding the risk-gene and mappability vetoes.
- 16/27 (59%) of Ehsan's candidates pass all of our hard vetoes.
- 11/27 would be excluded by our data — none of which a human-regulatory-element liftover can detect, since they depend on real canine population/structural data or canine gene identity:

| chrom | start | end | reason |
|---|---:|---:|---|
| NC_051812.1 | 5172690 | 5239184 | TAD boundary overlap |
| NC_051814.1 | 6647061 | 6715959 | TAD boundary overlap |
| NC_051816.1 | 729204 | 780179 | ATAC peak overlap |
| NC_051820.1 | 1623305 | 1690052 | ATAC peak overlap |
| NC_051828.1 | 16606256 | 16661497 | TAD boundary overlap |
| NC_051830.1 | 3084750 | 3142913 | TAD boundary overlap |
| NC_051843.1 | 6868917 | 6934242 | TAD boundary overlap |
| NC_051809.1 | 2620428 | 2685222 | risk gene nearby (OPCML — tumor suppressor) |
| NC_051812.1 | 6167980 | 6227889 | risk gene nearby (NOVA1) |
| NC_051823.1 | 15192889 | 15259450 | risk gene nearby (FAT4 — tumor suppressor) |
| NC_051843.1 | 8354252 | 8409529 | self-mappability check flagged (ambiguous realignment) |

## Update 2026-08-20 — full V2 cross-check (all 8 hard vetoes, not just V1's 4)

Re-checked all 27 against the current V2 pipeline (`candidates_scored_v5.tsv`
— miRNA, risk-gene radius, gene-dense neighborhood, lncRNA/smallRNA, TAD-content,
on top of the original TAD-boundary/ATAC-peak/risk-gene/mappability vetoes).

**Only 4/27 (15%) now survive**, down from 16/27 under the old V1-era vetoes:

| chrom | start | end | final_score |
|---|---:|---:|---:|
| NC_051807.1 | 10779807 | 10837319 | 0.642 |
| NC_051812.1 | 5534832 | 5600866 | 0.686 |
| NC_051821.1 | 4818014 | 4875945 | 0.487 |
| NC_051843.1 | 10578732 | 10643010 | 0.566 |

The former #1 (`NC_051812.1:52431-118675`) is still one of Ehsan's 27, but no
longer survives our V2: `veto_risk_gene_radius` (NLRP3 at 214kb) and
`veto_gene_dense_neighborhood`. Of the 23 now-excluded, every single one is
caught by `veto_gene_dense_neighborhood` (the third-gene-within-50kb check,
the strictest criterion in V2 and the one with no equivalent in Ehsan's
liftover-based filtering) either alone or combined with TAD boundary, ATAC
peak, risk-gene-radius, or TAD-content vetoes.

**Also checked Ehsan's newest, further-refined final 7** (his own second-round
filtering: 27→18→7, using real ATAC-seq from 5 canine PBMC samples + CanFam3.1
regulatory elements/CpG islands liftover, shared 2026-08-20): **0/7 survive
V2**, all 7 via `veto_gene_dense_neighborhood`.

## Update 2026-08-21 — V9 check against Ehsan's own regulatory-element set

Separately from candidate-list agreement above, Ehsan sent the actual
regulatory-element BED behind his own filtering (`ehsan_regulatory_elements_ROS.bed`,
75,600 elements, CanFam3.1 lifted to ROS_Cfam_1.0 by him) — added as a
formal hard veto in V9 (`score_ship_candidates.py`, `veto_external_regulatory_element`).

Of the 3 candidates reported as still standing on both sides (the most
recent recheck against V8, referenced in the 2026-08-21 email update — not
separately tabulated in this file before now): **2 of the 3 are directly
excluded by Ehsan's own regulatory-element set**:

| chrom | start | end | V9 status |
|---|---:|---:|---|
| NC_051807.1 | 10779807 | 10837318 | survives — rank 7 in V9, score 0.7369 |
| NC_051821.1 | 4818014 | 4875944 | **excluded** — overlaps Ehsan's regulatory elements (only veto tripped) |
| NC_051843.1 | 10578732 | 10643009 | **excluded** — overlaps Ehsan's regulatory elements (only veto tripped) |

Worth being precise about what this does and doesn't mean: candidate-list
convergence (both pipelines independently generating/keeping the same
locus) and regulatory-element overlap (a specific annotation-based check)
are different questions. A candidate can be a convergence hit and still
sit on a real regulatory element neither pipeline's *candidate-generation*
step was checking for. This isn't a contradiction — it's exactly the kind
of thing a hard veto based on the actual annotation file, rather than
final coordinate-list comparison alone, is supposed to catch.

## Update 2026-09-09 — V12 check against Ehsan's exact final 7 (email 2026-08-20)

Ehsan's email (2026-08-20) gives the exact coordinates of his final 7
pGSH candidates (27->18->7: ATAC reproducibility in >=4/5 samples, then
CanFam3.1 regulatory-element/CpG-island filtering, coords converted to
ROS_Cfam_1.0) - saved verbatim to `05_SHIP/ehsan_final7_pGSH.tsv`. He also
separately sent `dog_40_chromosomes_result_SHIP.txt` (in
`dog_40_chromosomes_result_SHIP (1).zip`) - checked in full 2026-09-09: this
is confirmed to be his own copy of the same raw SHIP output this
repository's `ship_raw_candidates.tsv` is parsed from (461 records, same
structure), **not** an intermediate 27 or 18 list - no filtering-stage
markers exist in that file. His full 27-candidate coordinate list was never
saved as a file in this repo, only referenced in earlier emails; only the
final 7 and the specific loci already itemized elsewhere in this document
are available to check against.

All 7 of his final candidates match **exactly** (identical coordinates) to
entries in our own 461-candidate set, confirming both pipelines share the
same SHIP run. All 7 are hard-vetoed by our full pipeline (V12,
2026-09-09), 6 of 7 via `veto_gene_dense_neighborhood` alone or combined
with others - consistent with the V2-era finding above (his refined-7 vs.
our `veto_gene_dense_neighborhood`), now reconfirmed after the 2026-09-09
ATAC methodology fix (gate+gain, peak_frequency-based scoring,
accessibility floor) changed nothing about this specific outcome, since
none of that touches the gene-density check:

| Ehsan's locus | Our exclusion reason(s) |
|---|---|
| NC_051807.1:73210443-73261055 | risk_gene_radius, gene_dense_neighborhood |
| NC_051812.1:788172-845815 | gene_dense_neighborhood |
| NC_051814.1:6647061-6715959 | tad_boundary, gene_dense_neighborhood |
| NC_051814.1:16586697-16654416 | gene_dense_neighborhood |
| NC_051823.1:22718008-22788982 | tad_boundary, atac_peak, gene_dense_neighborhood |
| NC_051843.1:50115742-50188659 | gene_dense_neighborhood, tad_risk_gene |
| NC_051843.1:56733328-56786162 | risk_gene_radius, gene_dense_neighborhood, tad_risk_gene |

The durable convergence point remains `NC_051807.1:10779807-10837318`
(from his *original* 27, not the refined 7) - survives every version V2
through V12 of our pipeline.

## Takeaway

The two independent approaches substantially agree at the V1 level (TAD/ATAC/
risk-gene/mappability only), including at the top of the ranking — a
meaningful in-silico validation signal in itself (Vasco's 2026-04-17 email
asked specifically for this kind of cross-check). Once the full Ahmed et al.
2026 checklist (V2) is applied, agreement drops sharply (27→4, 7→0), driven
almost entirely by one criterion — a third gene within 50kb of either window
edge — that has no equivalent in a human-regulatory-liftover approach, since
it depends on canine gene density/identity directly rather than homology.
This is now the central open discussion point with Ehsan (see the email
thread) rather than a solved cross-validation.
