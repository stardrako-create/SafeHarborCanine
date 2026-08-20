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
