# 04_tracks_processadas

Genome-wide tracks produced by the pipelines in `../scripts/`, one
subdirectory per assembly (`ROS_Cfam_1.0/`) and layer (`ATAC/`, `RRBS/`,
`HiC/`).

Small tabular summaries (QC weight tables, consensus peak calls) are
tracked directly in this repository. Large binary tracks — bigWig
(`.bw`), Hi-C contact matrices (`.mcool`), and the redundant bedgraph
export of the ATAC peak-frequency track — are excluded from git (several
exceed GitHub's 100 MB per-file limit) and are released as a versioned
dataset on Zenodo: **[10.5281/zenodo.22003934](https://doi.org/10.5281/zenodo.22003934)**.
Filenames there are flattened and prefixed by layer (`ATAC_`, `RRBS_`,
`HiC_`) since Zenodo has no subdirectories; the table below maps each
Zenodo filename to its local pipeline path and description.

| Zenodo filename | Local pipeline path | Layer | Description |
|---|---|---|---|
| — | `ATAC/mother_track/qc_weights.tsv` | ATAC | Per-dog QC confidence weight used in the Mother Track |
| — | `ATAC/mother_track/consensus_peaks_ATAC.bed` | ATAC | Consensus accessible-region peak calls across all 71 dogs |
| `ATAC_weighted_mean.bw` | `ATAC/mother_track/mother_track_weighted_mean.bw` | ATAC | Weighted-mean accessibility signal, genome-wide |
| `ATAC_peak_frequency.bw` / `.bedgraph` | `ATAC/mother_track/peak_frequency.bw` / `.bedgraph` | ATAC | Fraction of dogs with a called peak at each position |
| `ATAC_variability.bw` | `ATAC/mother_track/variability.bw` | ATAC | Signal variability across dogs |
| — | `RRBS/rrbs_qc_weights.tsv` | RRBS | Per-dog QC confidence weight used in the methylation Mother Track |
| `RRBS_weighted_mean.bw` | `RRBS/methylation_weighted_mean.bw` | RRBS | Weighted-mean CpG methylation (%), genome-wide |
| `RRBS_cpg_coverage_frequency.bw` | `RRBS/cpg_coverage_frequency.bw` | RRBS | CpG coverage frequency across dogs |
| `RRBS_variability.bw` | `RRBS/variability.bw` | RRBS | Methylation variability across dogs |
| `HiC_contact_matrix.mcool` | `HiC/mischka_hic.mcool` | Hi-C | Multi-resolution contact matrix (Mischka, 3 merged libraries) |
| `HiC_tad_insulation.tsv` | `HiC/tad_insulation.tsv` | Hi-C | Raw insulation score table (100/250/500kb windows) |
| `HiC_tad_insulation_100000.bw` / `_250000.bw` / `_500000.bw` | `HiC/tad_insulation.tsv.<window>.bw` | Hi-C | Genome-wide TAD insulation score bigwigs |
| — | `HiC/tad_boundaries.bed` | Hi-C | Called TAD boundary positions (8,367 calls) |

Rows marked "—" are small enough to be tracked directly in git and don't
have a separate Zenodo copy.
