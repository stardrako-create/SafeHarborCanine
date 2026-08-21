# Consensus ATAC peak threshold — validation against independent regulatory annotation

Ehsan asked (2026-08-20 email) whether the consensus-peak veto (≥36/71 dogs,
ArchR-style majority threshold) actually identifies functional regulatory
regions, rather than being trusted by convention alone. Tested empirically
against two independent, dog-native regulatory proxies — chosen to avoid
relying on a cross-species liftover of human/mouse regulatory annotations,
which is unreliable for enhancers specifically at this evolutionary distance:

- **TSS ± 2kb** of every annotated gene in the ROS_Cfam_1.0 GFF3 (37,061
  regions, strand-aware, derived directly from `genomic.gff` — a standard
  promoter-proximal proxy).
- **CpG islands**, called directly on the ROS_Cfam_1.0 sequence itself
  (Gardiner-Garden & Frommer 1987 criteria: ≥200bp window, GC% ≥50%,
  Obs/Exp CpG ≥0.6 — the same definition behind UCSC's `cpgIslandExt`
  track), 284,435 islands genome-wide. Computed natively rather than
  lifted from another species, so it carries no cross-species assumption.

## Method

`bedtools intersect` between `consensus_peaks_ATAC.bed` (8,840 peaks) and
each proxy, compared against each proxy's genome-wide base-pair coverage
fraction as the expected background rate under no enrichment.

## Result

| Proxy | Peaks overlapping | Observed fraction | Genome background | Enrichment |
|---|---:|---:|---:|---:|
| TSS ±2kb | 5,941 / 8,840 | 67.2% | 6.18% | **10.9x** |
| CpG islands | 6,656 / 8,840 | 75.3% | 4.98% | **15.1x** |
| Either | 6,962 / 8,840 | 78.8% | 11.17% (upper bound) | **7.0x** (conservative) |

## Conclusion

The ≥36/71-dog consensus threshold is not an arbitrary statistical
convention borrowed from ArchR without local justification — it identifies
regions 7-15x enriched for independent, dog-native evidence of regulatory
function (promoter proximity and CpG island content) relative to genome-wide
background. This is real support for treating a consensus ATAC peak as
"looks like an actual regulatory element," which is the entire justification
for using it as a hard veto rather than folding it into the soft score.

No change made to the threshold itself — this was a validation check, not a
recalibration; the enrichment supports keeping ≥36/71 as-is.
