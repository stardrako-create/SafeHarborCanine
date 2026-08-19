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

## Result

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

## Takeaway

The two independent approaches substantially agree, and agree at the top of
the ranking — a meaningful in-silico validation signal in itself (Vasco's
2026-04-17 email asked specifically for this kind of cross-check). Where
they disagree, it is because this repository's filter is grounded in real
canine population/structural data and canine gene identity rather than
human homology. The 11 divergent candidates are a concrete, explainable set
to discuss before finalizing a shared shortlist — the addition of the
risk-gene and mappability checks (2026-08-19) caught 4 more (OPCML, NOVA1,
FAT4 proximity; one ambiguous self-alignment) that the original TAD/ATAC-only
version of this pipeline had missed.
