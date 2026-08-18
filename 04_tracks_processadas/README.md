# 04_tracks_processadas

Genome-wide tracks produced by the pipelines in `../scripts/`, one
subdirectory per assembly (`ROS_Cfam_1.0/`) and layer (`ATAC/`, `RRBS/`,
`HiC/`).

Small tabular summaries (QC weight tables, consensus peak calls) are
tracked directly in this repository. Large binary tracks — bigWig
(`.bw`), Hi-C contact matrices (`.mcool`), and the redundant bedgraph
export of the ATAC peak-frequency track — are excluded from git (several
exceed GitHub's 100 MB per-file limit) and are released as a versioned
dataset on Zenodo. The DOI will be added here once published.

| File | Layer | Description |
|---|---|---|
| `ATAC/mother_track/qc_weights.tsv` | ATAC | Per-dog QC confidence weight used in the Mother Track |
| `ATAC/mother_track/consensus_peaks_ATAC.bed` | ATAC | Consensus accessible-region peak calls across all 71 dogs |
| `ATAC/mother_track/mother_track_weighted_mean.bw` * | ATAC | Weighted-mean accessibility signal, genome-wide |
| `ATAC/mother_track/peak_frequency.bw` * | ATAC | Fraction of dogs with a called peak at each position |
| `ATAC/mother_track/variability.bw` * | ATAC | Signal variability across dogs |
| `RRBS/rrbs_qc_weights.tsv` | RRBS | Per-dog QC confidence weight used in the methylation Mother Track |
| `RRBS/methylation_weighted_mean.bw` * | RRBS | Weighted-mean CpG methylation (%), genome-wide |
| `RRBS/cpg_coverage_frequency.bw` * | RRBS | CpG coverage frequency across dogs |
| `RRBS/variability.bw` * | RRBS | Methylation variability across dogs |
| `HiC/*.mcool` * | Hi-C | Multi-resolution contact matrix (Mischka, 3 merged libraries) |
| `HiC/tad_insulation_score.bw` * | Hi-C | Genome-wide TAD insulation score |
| `HiC/tad_boundaries.bed` | Hi-C | Called TAD boundary positions |

\* not tracked in git — see Zenodo release.
