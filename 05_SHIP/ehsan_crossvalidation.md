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
population data (71-dog ATAC/RRBS, single-individual Hi-C) rather than
human-liftover homology, narrowing the same 461 → 280 ranked candidates.

## Result

- All 27 of Ehsan's candidates are confirmed present in the shared 461-candidate set.
- **Our #1-ranked candidate (`NC_051812.1:52431-118675`, score 0.83) is also in Ehsan's final 27** — two independent filtering strategies converge on the same top locus.
- 20/27 (74%) of Ehsan's candidates pass our hard vetoes too.
- 7/27 would be excluded by our data: 5 overlap a Hi-C TAD boundary, 2 overlap
  a real canine ATAC consensus peak. Neither of these can be detected by a
  human-regulatory-element liftover, since they reflect canine-specific
  population epigenomic/structural signal, not homology to human annotation.

| chrom | start | end | our verdict | reason |
|---|---:|---:|---|---|
| NC_051812.1 | 5172690 | 5239184 | excluded | TAD boundary overlap |
| NC_051814.1 | 6647061 | 6715959 | excluded | TAD boundary overlap |
| NC_051816.1 | 729204 | 780179 | excluded | ATAC peak overlap |
| NC_051820.1 | 1623305 | 1690052 | excluded | ATAC peak overlap |
| NC_051828.1 | 16606256 | 16661497 | excluded | TAD boundary overlap |
| NC_051830.1 | 3084750 | 3142913 | excluded | TAD boundary overlap |
| NC_051843.1 | 6868917 | 6934242 | excluded | TAD boundary overlap |

## Takeaway

The two independent approaches substantially agree, and agree at the top of
the ranking — a meaningful in-silico validation signal in itself (Vasco's
2026-04-17 email asked specifically for this kind of cross-check). Where
they disagree, it is because this repository's filter is grounded in real
canine data rather than human homology, which is the case for a dog-specific
project. The 7 divergent candidates are a concrete, explainable set to
discuss before finalizing a shared shortlist.
