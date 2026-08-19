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

## Re-evaluation against V2 (Ahmed et al. 2026 checklist criteria, 2026-08-19)

| Check | Result |
|---|---|
| miRNA within 300kb? | No — nearest is 31 Mb away. Passes. |
| Risk gene within 300kb radius? | **Fails** — `FRK` (a risk gene) sits 129 kb away |
| Overlaps lncRNA/small RNA? | No — passes |
| Outside transcriptional unit? | **Fails** — the window is fully contained inside `NT5DC1` (`NC_051816.1:72,672,269-72,791,171`) on ROS_Cfam_1.0. This locus is intragenic here, not intergenic — a genuinely different situation from the SHIP candidates, which are intergenic by construction. (Precedent exists for intragenic safe harbors — AAVS1 sits in PPP1R12C, hROSA26 in THUMPD3 — so this alone isn't disqualifying, but it does mean criterion 4 fails outright.) |
| Own TAD risk-gene content | Could not resolve — this point doesn't fall inside any derived TAD interval (likely sits within a boundary block itself) |

## Conclusion

Combined with the already-elevated ATAC accessibility (0.214, vs. 0.04–0.11
typical for passing SHIP candidates) noted below, this candidate **would
not pass V2**: it fails the risk-gene-radius check (FRK, 129 kb) and is
intragenic in `NT5DC1`, not intergenic. It was found by a different,
non-systematic method (manual JBrowse inspection on a different assembly),
so it remains a useful independent data point and a reminder of why
systematic scoring matters — but it should not be treated as equivalent to
the top SHIP-derived V2 candidate (`NC_051811.1:48,020,921-48,077,046`).

## V1-era evaluation (superseded by the V2 checks above, kept for context)

| Check | Result |
|---|---|
| One of the 461 SHIP candidates? | **No** — nearest SHIP candidate on this chromosome is ~370 kb away. Found by manual inspection, not SHIP's intergenic-interval search, so non-membership isn't itself disqualifying. |
| Overlaps a Hi-C TAD boundary? | No — passes |
| Overlaps an ATAC consensus peak? | No — passes |
| ATAC weighted-mean accessibility | 0.214 — notably higher than passing SHIP candidates (typically 0.04–0.11). Below our peak-calling threshold, but flagged: elevated accessibility this far above the candidate population could indicate proximity to an unannotated regulatory element, even without a called peak. |
| RRBS weighted-mean methylation | 2.4% — favorable, comparable to the better-ranked SHIP candidates |
| Distance to nearest TAD boundary | ~93 kb — reasonable, well inside a TAD |
