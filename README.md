# Safe Harbor Canino

[![DOI (code)](https://zenodo.org/badge/DOI/10.5281/zenodo.21996452.svg)](https://doi.org/10.5281/zenodo.21996452)
[![DOI (data)](https://zenodo.org/badge/DOI/10.5281/zenodo.22003933.svg)](https://doi.org/10.5281/zenodo.22003933)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Identification of genomic safe harbor loci in the domestic dog (*Canis lupus
familiaris*), integrating population-level chromatin accessibility, DNA
methylation, and single-individual 3D genome structure, for use as candidate
integration sites in canine CAR-T cell engineering.

Developed in the lab of Vasco M. Barreto (NOVA FCT). See
[CITATION.cff](CITATION.cff) for the authors.

## Status at a glance — 12 October 2026

The defined computational review is complete for the data currently available.
The shortlist is **w01** plus one additional provisional candidate,
**ANO2-39708**. These are computational candidates: none has been validated
experimentally for integration, expression or function, and no guide RNAs are
claimed.

| Candidate | ROS_Cfam_1.0 interval (BED, 0-based, end-exclusive) | Disposition |
|---|---|---|
| w01 | NC_051811.1:17255706–17256706 | First priority for testing |
| ANO2-39708 | NC_051831.1:39708074–39709074 | Additional provisional candidate |
| w11 | NC_051836.1:27041161–27042161 | Not advanced: structural context unresolved |
| bg1k_2848, bg1k_0442 | NC_051840.1:28808728–28809728, NC_051837.1:10430486–10431486 | Genomic comparators, not validated negative controls |

Start with [the snapshot guide](RELEASE_20261012.md), then
[the full results and inventories](08_CANDIDATE_CLOSURE_AUDIT/FINAL_RESULTS_20261012.md)
and [the decision log](08_CANDIDATE_CLOSURE_AUDIT/CURRENT_STATE_AND_DECISIONS.md)
(working notes are in Portuguese). Earlier snapshots:
[25 September 2026](RELEASE_20260925.md). Older version directories preserve the
development history; their candidate rankings are historical.

### Main limitations

- Everything here is computational. Integration, expression and function have not been tested, and no guide RNAs or nuclease are specified.
- The population layers are PBMC from a mixed-cell blood cohort, not purified T cells, and direct T-cell accessibility support at the shortlisted windows is sparse.
- Methylation coverage is incomplete (RRBS observes 20% of the w01 window and 15% of the ANO2-39708 window), and the risk-gene annotation is incomplete (53 of 2,654 symbols unresolved for ANO2-39708).
- w11 was withdrawn because its structural context could not be resolved with short reads; this does not show that the locus is unsafe.
- A complete end-to-end rerun of the later audit stages has not been performed.

## Background

A genomic safe harbor is a locus where a transgene can be integrated with a
low risk of disrupting endogenous gene function or triggering oncogenic
activation — a prerequisite for safely engineering canine CAR-T cells.
Identifying such loci requires converging evidence across several independent
layers of genome biology: which regions are open and transcriptionally
permissive across individuals, which regions carry a stable epigenetic
signature, and how the genome folds in 3D (so a candidate site isn't sitting
on a topologically associating domain (TAD) boundary or chromatin loop
anchor whose disruption could deregulate a distant gene).

This repository holds the full analysis pipeline — not a single script, but
the reproducible chain from public raw sequencing data to genome-wide
annotation tracks — used to build that converging evidence for the dog
reference genome ROS_Cfam_1.0 (GCF_014441545.1).

## Data layers

| Layer | Assay | Cohort | Source | Status |
|---|---|---|---|---|
| Chromatin accessibility | ATAC-seq | 71 dogs, PBMC | Jin et al. 2024, *Aging Cell* — [PRJNA1048909](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1048909) | Done |
| DNA methylation | RRBS | 71 dogs, PBMC | Jin et al. 2024, *Aging Cell* — [PRJNA1049514](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1049514) | Done |
| 3D genome structure | Hi-C | 1 dog ("Mischka"), 3 assembly Hi-C libraries; cell type not reported | Wang et al. 2021, *Communications Biology* — [PRJNA587469](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA587469) | Done |

The ATAC and RRBS layers are population-level: 71 individuals from the same
PBMC cohort, each contributing a QC-weighted vote to a "Mother Track" that
summarizes accessibility / methylation genome-wide with a confidence measure
at every base pair. The Hi-C layer is structural context from a single dog
of a different breed; the paper reports blood-derived HMW DNA but does not
explicitly report the cell/nuclear preparation used for the Hi-C libraries — its
Dovetail Hi-C library is what originally scaffolded the **UU_Cfam_GSD_1.0**
assembly, which ROS_Cfam_1.0 itself was chromosome-scaffolded against via
RagTag, giving it strong chromosome-level synteny with the reference used
throughout this project. Because it is a single non-cohort individual with
an imperfect tissue match, **the Hi-C layer is intentionally weighted low
relative to the ATAC/RRBS Mother Tracks** in any downstream safe-harbor
scoring — it is used as structural context (avoid known TAD boundaries /
loop anchors), not as a population-level confidence signal.

Later versions also use a 76-dog ATAC cohort (the 71 dogs plus 5 additional
healthy control dogs) built in a single unified run.

## Repository structure

```
01_referencia/                   Reference genome + annotation (ROS_Cfam_1.0). Not tracked — public, see below.
02_dados_brutos/                 Raw SRA/FASTQ. Not tracked — regenerable from the accessions in scripts/config*.yaml.
03_tracks_originais/             Third-party tracks used as-is (e.g. CanFam3.1 liftover sources).
04_tracks_processadas/           Small tables and peak sets; large bigWig/.mcool tracks are released via Zenodo (see below).
05_SHIP/                         SHIP candidate search and per-locus scoring (V1–V12), including the Ahmed et al. 2026
                                 checklist comparison (ahmed2026_checklist_comparison.md) and VERSIONS.md.
06_v1.20.0 … 06_v1.23.1/         Version-by-version rebuild and audit stages; most have a CURRENT.md (Portuguese).
07_FINAL_CANDIDATES_2026-09-19/  Candidate panel of 19 September 2026 (superseded by 08_).
08_CANDIDATE_CLOSURE_AUDIT/      Local-window audit, structural-variant review, ANO2-39708 review and the
                                 final results of 12 October 2026 (FINAL_RESULTS_20261012.md).
scripts/                         All pipeline code: Snakemake workflows, per-layer processing scripts, configs.
CURRENT.md                       Chronological working log (Portuguese).
06_GEG-SH/, 07_integracao/, 08_candidatos_finais/
                                 Placeholders from the initial layout; the work moved to the folders above.
```

Each data directory that isn't tracked in git has its own `README.md`
explaining what belongs there and how to regenerate it.

## Pipelines

All three layers follow the same pattern: align raw reads → per-individual
(or per-library) QC and filtering → per-individual confidence weighting →
merge into a genome-wide track. Full detail is in the scripts themselves;
summary below.

### ATAC-seq (`scripts/Snakefile_persample.smk`, `config.yaml`)
Bowtie2 alignment → mitochondrial/blacklist filtering and deduplication →
MACS3 peak calling per dog → per-dog QC scoring (`compute_qc_weights.py`) →
weighted Mother Track and consensus peak set. The unified builder
`build_mother_track_v2.py` (all dogs in one run, confidence gate plus
fold-enrichment units) supersedes `build_mother_track.py` and
`join_atac_cohorts.py`, which are kept for history; `call_consensus_peaks.py`
builds the consensus peaks. The primary peak-discriminative signal is the
per-bin *peak frequency* (how many dogs have their own independent MACS3 peak
nearby); the weighted-mean track is used only as a coarse accessibility level.

### RRBS (`scripts/Snakefile_rrbs.smk`, `config_rrbs.yaml`)
Trim Galore → Bismark bisulfite alignment and per-CpG methylation extraction
→ per-dog QC scoring (`qc_metrics_rrbs.py`) → weighted methylation Mother
Track (mean methylation, coverage frequency, variability;
`build_methylation_track.py`).

### Hi-C (`scripts/Snakefile_hic.smk`, `config_hic.yaml`)
`bwa mem -5SP -T0` chimeric-aware alignment (4DN/distiller-nf convention) →
`pairtools parse` → `sort` → `dedup` → `select` (`UU`, MAPQ ≥ 30; no
restriction-fragment filtering, since the Dovetail kit's enzyme was never
disclosed in the source publication) per library → merge across the 3
libraries of the same individual → `cooler cload`/`zoomify`/`balance` →
`cooltools insulation` for TAD boundary calls (`build_hic_tracks.py`).

### Candidate scoring and audit (`scripts/score_ship_candidates_v2.py`, `05_SHIP/`, `08_CANDIDATE_CLOSURE_AUDIT/`)
The 461 candidate regions from the SHIP search are scored against hard vetoes
derived from the eight-criterion checklist of Ahmed et al. 2026 (plus
project-specific additions) and a continuous score; thresholds taken from the
literature (for example the 50-kb gene-distance criterion) are kept at their
published values. Regions that survive are then audited as 1-kb local windows
(repeats, structural-variant catalogues, conservation, variants, risk-gene
annotation), with every evidence gap recorded rather than treated as a pass. Per-version changes are documented in
`05_SHIP/VERSIONS.md` and the `CURRENT.md` file of each version folder.

Every per-dog/per-library QC weight, and the exact formula combining local
sequencing depth confidence with global sample-quality metrics, is computed
transparently in the corresponding `compute_qc_weights.py` / QC scripts —
nothing about "which individuals count more" is a black box.

### Reproducing a run

```bash
conda env create -f environment.yml
conda activate atac
cd scripts

# ATAC-seq
snakemake -s Snakefile_persample.smk --configfile config.yaml --cores <N>

# RRBS
snakemake -s Snakefile_rrbs.smk --configfile config_rrbs.yaml --cores <N>

# Hi-C
snakemake -s Snakefile_hic.smk --configfile config_hic.yaml --cores <N> --resources wip=1
python build_hic_tracks.py --config config_hic.yaml
```

Sample accessions and all local paths are declared in the `config*.yaml`
files — adjust `paths:` to your own filesystem before running. A complete
end-to-end rerun of the later audit stages has not been performed; see the
limits stated in [RELEASE_20261012.md](RELEASE_20261012.md).

## Data availability

Raw sequencing data is not redistributed here: it is already public under
the SRA BioProject accessions listed above, and the exact sample list used
is in `scripts/config.yaml` / `config_rrbs.yaml` / `config_hic.yaml`. The
reference genome (ROS_Cfam_1.0, GCF_014441545.1) is public via NCBI.

Processed genome-wide tracks (bigWig coverage/methylation/insulation tracks,
`.mcool` Hi-C contact matrices) are too large for git (several files exceed
GitHub's 100 MB per-file limit) and are released as a versioned dataset on
Zenodo: **[10.5281/zenodo.22003933](https://doi.org/10.5281/zenodo.22003933)**
(concept DOI, always resolves to the latest version). The archived tracks come
from the original per-cohort pipeline (71-dog, 5-dog and joined 76-dog ATAC
tracks); the later unified rebuild is documented in the code, but its large
output tracks are not archived. Small tabular outputs (QC weight tables, candidate loci
lists) are tracked directly in this repository under `04_tracks_processadas/`
and `05_SHIP/` through `08_CANDIDATE_CLOSURE_AUDIT/`.

## Related work

[LocusAtlas](https://github.com/stardrako-create/EpiLog) (formerly EpiLog; [live catalog](https://locusatlas.netlify.app))
is a literature-derived catalog of known genomic safe harbor loci across
species, built by the same author. Its dog-species survey found essentially
no published safe-harbor characterization for *Canis lupus familiaris* —
that gap is what this repository addresses directly, with original
population-level ATAC-seq/RRBS data and single-individual Hi-C structural
context, rather than literature mining.

## Contributions and credits

- **Manuel Paiva Sequeira**: design, implementation and execution of the pipeline (ATAC-seq, RRBS and Hi-C processing, the Mother Track construction, candidate scoring, the local-window audit and the structural-variant, conservation and variant reviews), data curation, documentation and review of the results.
- **Ehsan Valiollahi**: the initial idea and the initial SHIP candidate generation (the 461 convergent intergenic regions on ROS_Cfam_1.0). He also provided an independent regulatory-element set and five additional control ATAC-seq samples used in later versions, and his review caught two errors that are documented in the dataset record.
- **Vasco M. Barreto**: supervision; principal investigator of the laboratory in which the project is carried out (NOVA FCT).

Third-party methods and data (SHIP, GEG-SH, the checklist of Ahmed et al. 2026 and the public ATAC-seq, RRBS and Hi-C data sets) are credited under Citation below.

**AI assistance.** Parts of the code, analysis and documentation were produced with AI coding assistants (Anthropic's Claude Code and OpenAI's Codex) working under the author's direction. The scripts, receipts, regression tests and the SHA-256 manifest in this repository (`python scripts/verify_manifest.py PUBLICATION_MANIFEST_20261012.json`) are the verifiable record of what was run.

## License

MIT License — see [LICENSE](LICENSE).

## Citation

This repository is archived on Zenodo; the concept DOI
[10.5281/zenodo.21996452](https://doi.org/10.5281/zenodo.21996452) always
resolves to the latest archived snapshot of the code, and
[10.5281/zenodo.22003933](https://doi.org/10.5281/zenodo.22003933) to the
latest version of the dataset. Machine-readable citation metadata is in
[CITATION.cff](CITATION.cff) (use GitHub's "Cite this repository" button, top
right of the repo page). A dedicated preprint reference will be added once
available. In the meantime, please also cite the source datasets where
relevant:

- Jin et al. 2024, *Aging Cell* (ATAC-seq / RRBS PBMC cohort, PRJNA1048909 / PRJNA1049514)
- Wang et al. 2021, *Communications Biology* (Hi-C / UU_Cfam_GSD_1.0 assembly, PRJNA587469)
- Leitão et al. 2025, *Scientific Reports*, 15, 7193 — [SHIP](https://github.com/MCLeitao/Ship), the candidate-generation tool `05_SHIP/` is built on (DOI: [10.1038/s41598-025-91249-9](https://doi.org/10.1038/s41598-025-91249-9))
- Shrestha et al. 2022, *Genome Biology*, 23, 199 — [GEG-SH](https://github.com/dewshr/GEG-SH), whose tissue-specific hard-veto/soft-score framework `05_SHIP/score_ship_candidates.py` follows (DOI: [10.1186/s13059-022-02770-3](https://doi.org/10.1186/s13059-022-02770-3))
- Ahmed et al. 2026, *Cells*, 15, 81 — comprehensive review of human SHS selection criteria; `05_SHIP/ahmed2026_checklist_comparison.md` scores this pipeline against its 8-criterion checklist (DOI: [10.3390/cells15010081](https://doi.org/10.3390/cells15010081))
