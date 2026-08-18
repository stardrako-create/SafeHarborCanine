# The original chr12 candidate, lifted to ROS_Cfam_1.0

The project's first candidate (visual JBrowse/iDog scouting, documented in
"Identificação de um Possível Safe Harbor Canino para Aplicação em CAR.pdf")
was identified at `chr12:72,350,025-72,351,060` on **UU_Cfam_GSD_1.0**
(`NC_049233.1`, iDog assembly parameter confirms this) — not ROS_Cfam_1.0,
the assembly this repository's pipeline and the SHIP candidate lists use.
This mapping was proposed but never completed in the Manuel↔Ehsan email
thread (2026-04-29): "I would therefore first map my chr12 sequence to the
same reference and then check whether it passes or fails each SHIP filter."

## Liftover method

No chain file exists between these two independently-assembled dog
genomes, so the region was lifted by direct realignment:
1. Fetched the exact 1,036 bp sequence at `NC_049233.1:72,350,025-72,351,060`
   via NCBI E-utilities (`efetch`).
2. Aligned it with `bwa mem` against the ROS_Cfam_1.0 genome (the same
   index built for the Hi-C pipeline).

Result: a near-perfect match (MAPQ 60, 1036M CIGAR, only 2 mismatches out
of 1036 bp) at **`NC_051816.1:72,767,426-72,768,461`**.

## Evaluation against this repository's data

| Check | Result |
|---|---|
| One of the 461 SHIP candidates? | **No** — nearest SHIP candidate on this chromosome is ~370 kb away. Found by manual inspection, not SHIP's intergenic-interval search, so non-membership isn't itself disqualifying. |
| Overlaps a Hi-C TAD boundary? | No — passes |
| Overlaps an ATAC consensus peak? | No — passes |
| ATAC weighted-mean accessibility | 0.214 — notably higher than passing SHIP candidates (typically 0.04–0.11). Below our peak-calling threshold, but flagged: elevated accessibility this far above the candidate population could indicate proximity to an unannotated regulatory element, even without a called peak. |
| RRBS weighted-mean methylation | 2.4% — favorable, comparable to the better-ranked SHIP candidates |
| Distance to nearest TAD boundary | ~93 kb — reasonable, well inside a TAD |

## Conclusion

The original chr12 candidate passes both of this repository's hard vetoes
and has favorable methylation and TAD-boundary distance, but its
accessibility is elevated relative to the SHIP-derived candidate
population — worth a closer look (e.g. inspect the ATAC signal shape
directly rather than just the interval mean) before treating it as
equally strong as the top SHIP-derived candidates in `candidates_scored.tsv`.
It was found by a different, non-systematic method, so it is a useful
independent data point rather than a member of the ranked list.
