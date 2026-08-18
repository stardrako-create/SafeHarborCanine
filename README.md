# Safe Harbor Canino

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21996453.svg)](https://doi.org/10.5281/zenodo.21996453)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Identification of genomic safe harbor loci in the domestic dog (*Canis lupus
familiaris*), integrating population-level chromatin accessibility, DNA
methylation, and single-individual 3D genome structure, for use as candidate
integration sites in canine CAR-T cell engineering.

Maintained by the Vasco M. Barreto lab.

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
| 3D genome structure | Hi-C | 1 dog ("Mischka"), blood | Wang et al. 2021, *Communications Biology* — [PRJNA587469](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA587469) | In progress |

The ATAC and RRBS layers are population-level: 71 individuals from the same
PBMC cohort, each contributing a QC-weighted vote to a "Mother Track" that
summarizes accessibility / methylation genome-wide with a confidence measure
at every base pair. The Hi-C layer is structural context from a single dog
of a different breed and tissue (whole blood, not isolated PBMCs) — its
Dovetail Hi-C library is what originally scaffolded the **UU_Cfam_GSD_1.0**
assembly, which ROS_Cfam_1.0 itself was chromosome-scaffolded against via
RagTag, giving it strong chromosome-level synteny with the reference used
throughout this project. Because it is a single non-cohort individual with
an imperfect tissue match, **the Hi-C layer is intentionally weighted low
relative to the ATAC/RRBS Mother Tracks** in any downstream safe-harbor
scoring — it is used as structural context (avoid known TAD boundaries /
loop anchors), not as a population-level confidence signal.

## Repository structure

```
01_referencia/          Reference genome + annotation (ROS_Cfam_1.0). Not tracked — public, see below.
02_dados_brutos/        Raw SRA/FASTQ. Not tracked — regenerable from the accessions in scripts/config*.yaml.
03_tracks_originais/    Third-party tracks used as-is (e.g. CanFam3.1 liftover sources).
04_tracks_processadas/  Genome-wide bigWig/.mcool tracks produced by this pipeline. Not tracked in git — released via Zenodo/GitHub Releases (see below).
05_SHIP/                Safe Harbor Integration Prioritization: per-locus scoring that combines all layers.
06_GEG-SH/              Gene/element-context annotation of candidate loci.
07_integracao/          Cross-layer integration and candidate ranking.
08_candidatos_finais/   Final candidate safe harbor loci.
scripts/                All pipeline code: Snakemake workflows, per-layer processing scripts, configs.
logs/                   Local run logs. Not tracked.
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
MACS2-style peak calling per dog → per-dog QC scoring
(`compute_qc_weights.py`) → weighted Mother Track and consensus peak set
(`build_mother_track.py`, `call_consensus_peaks.py`).

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
files — adjust `paths:` to your own filesystem before running.

## Data availability

Raw sequencing data is not redistributed here: it is already public under
the SRA BioProject accessions listed above, and the exact sample list used
is in `scripts/config.yaml` / `config_rrbs.yaml` / `config_hic.yaml`. The
reference genome (ROS_Cfam_1.0, GCF_014441545.1) is public via NCBI.

Processed genome-wide tracks (bigWig coverage/methylation/insulation tracks,
`.mcool` Hi-C contact matrices) are too large for git (several files exceed
GitHub's 100 MB per-file limit) and are released as versioned datasets on
Zenodo, linked here once published. Small tabular outputs (QC weight
tables, candidate loci lists) are tracked directly in this repository under
`04_tracks_processadas/`, `05_SHIP/` through `08_candidatos_finais/`.

## Status

ATAC-seq and RRBS Mother Tracks are complete for all 71 dogs. The Hi-C
pipeline is running; two of three libraries are fully processed. Layers
05–08 (scoring, annotation, integration, final candidates) begin once the
Hi-C contact matrix and TAD boundary calls are complete.

## Related work

[EpiLog](https://github.com/stardrako-create/EpiLog) ([live catalog](https://epilogbio.netlify.app))
is a literature-derived catalog of known genomic safe harbor loci across
species, built by the same author. Its dog-species survey found essentially
no published safe-harbor characterization for *Canis lupus familiaris* —
that gap is what this repository addresses directly, with original
population-level ATAC-seq/RRBS data and single-individual Hi-C structural
context, rather than literature mining.

## License

MIT License — see [LICENSE](LICENSE).

## Citation

This repository is archived on Zenodo with a DOI covering the code as of
each tagged release: [10.5281/zenodo.21996453](https://doi.org/10.5281/zenodo.21996453).
Machine-readable citation metadata is in [CITATION.cff](CITATION.cff) (use
GitHub's "Cite this repository" button, top right of the repo page). A
dedicated preprint reference will be added once available. In the meantime,
please also cite the source datasets where relevant:

- Jin et al. 2024, *Aging Cell* (ATAC-seq / RRBS PBMC cohort, PRJNA1048909 / PRJNA1049514)
- Wang et al. 2021, *Communications Biology* (Hi-C / UU_Cfam_GSD_1.0 assembly, PRJNA587469)
