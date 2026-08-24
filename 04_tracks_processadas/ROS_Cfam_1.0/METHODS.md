# Methods — how each track in this dataset was actually built

This document exists so the Mother Track methodology is unambiguous, since
it is not a standard off-the-shelf pipeline — it's a specific weighting and
confidence scheme designed for this project. Every formula below is
transcribed directly from the scripts that produced the files (see the
named script for each section), not paraphrased from memory.

## ATAC-seq (71 dogs, PBMC, Jin et al. 2024, PRJNA1048909)

Scripts: `scripts/compute_qc_weights.py`, `scripts/build_mother_track.py`,
`scripts/call_consensus_peaks.py`.

**Per-dog global QC weight** (`compute_qc_weights.py`): for each dog *i*,
FRiP (fraction of reads in peaks) and TSS enrichment are each min-max
normalized to [0,1] across all 71 dogs, averaged, then rescaled onto
`[0.2, 1.0]` so no dog is ever fully zeroed out regardless of QC:

```
composite_i = (FRiP_norm_i + TSS_norm_i) / 2
weight_i    = 0.2 + 0.8 * composite_i
```

**Per-bin local confidence** (`build_mother_track.py`): a dog's raw
(non-normalized) coverage in a genomic bin is compared against that same
dog's own genome-wide mean raw coverage (λ_i, its background rate):

```
confidence_i(bin) = raw_i(bin) / (raw_i(bin) + λ_i)
```

This lets a low-depth dog's small raw count still register as real signal,
while an identical raw count from a high-depth dog (whose λ is much higher)
is recognized as noise-floor and down-weighted — and lets bins a given dog
simply wasn't measured at drop out of that dog's contribution instead of
being averaged in as "0% accessible."

**ATAC_weighted_mean.bw** (the main Mother Track), per bin:

```
weighted_mean(bin) = Σ_i [ weight_i · confidence_i(bin) · CPM_i(bin) ]
                      -------------------------------------------------
                      Σ_i [ weight_i · confidence_i(bin) ]
```

summed only over dogs with nonzero evidence at that bin.

**ATAC_variability.bw**: interquartile range (75th − 25th percentile) of
each dog's CPM at that bin, computed only over dogs whose local confidence
at that bin is ≥ 0.2 (the "confidence floor") — so a bin isn't scored as
highly variable just because most dogs weren't actually measured there.

**ATAC_peak_frequency.bw**: per-bin count of how many of the 71 dogs have
their own independently-called MACS3 narrowPeak overlapping that bin
(`bedtools multiinter` across all 71 per-dog peak sets). Context, not a
filter on its own.

**Consensus ATAC peaks** (`call_consensus_peaks.py`, used as the hard-veto
input in `05_SHIP/scripts/score_ship_candidates.py`): a region is called a
consensus peak if it is supported by **≥ ⌈(71+1)/2⌉ = 36 of the 71 dogs**
(strict majority — following ArchR's documented `addReproduciblePeakSet()`
convention, `reproducibility = "(n+1)/2"`), then regions within 75 bp of
each other are merged (`bedtools merge -d 75`). A percentile cutoff on the
continuous weighted-mean signal was tried first and rejected: that
distribution is narrow and non-bimodal (p50 = 0.09, p99.9 = 0.83) with no
natural peak/background separation — vote count across dogs does have one.

## ATAC-seq, extended cohort (76 dogs = 71 + 5, "ATAC_joined76") — 2026-08-24

Ehsan Valiollahi sent 5 additional canine ATAC-seq samples (GEO GSE278027,
SRR30799901–SRR30799905, PBMC, 2 breeds — Maltese ×4, Shih-Tzu ×1). GEO
metadata for all 5 traces to dogs "clinically diagnosed with mammary
tumors" (the study's growth protocol covers both healthy and tumor-bearing
animals; these 5 specifically are from the tumor-bearing arm) — a
different population from the healthy Jin et al. 2024 cohort, not a
straightforward addition of 5 more equivalent dogs. Processed through the
identical per-sample pipeline as the 71 (`Snakefile_persample.smk`: fastp
→ bowtie2 `--very-sensitive` → dedup/filter → MACS3 → QC metrics),
reusing the same ROS_Cfam_1.0 bowtie2 index.

**Standalone 5-dog Mother Track** (`ATAC_ehsan5/mother_track/`): same
`compute_qc_weights.py` / `build_mother_track.py` formulas above, applied
to just these 5 — min-max normalization for the QC weight is therefore
relative to this 5-dog cohort, not the 71-dog scale (weights: 0.200,
0.687, 1.000, 0.627, 0.326 for SRR30799901/02/03/04/05 respectively —
0.200 is the floor, not a computed near-zero value). Published as its own
track set, not merged blindly into the population baseline, given the
disease-context caveat above.

**Joined 76-dog track** (`ATAC_joined76/`): the two Mother Tracks combined
by **cohort size**, not by re-pooling all 76 raw dogs as equally-weighted
individuals:

```
joined(bin) = [ 71 · mother_mean_71(bin)  +  5 · mother_mean_5(bin) ] / 76
```

applied per-bin to both `mother_track_weighted_mean` and `variability`
(25 bp bins, matching the source tracks). This is the standard "pooled
mean = size-weighted average of group means" identity, applied explicitly
at the group level since each cohort's own QC-weighted average was already
computed independently — verified numerically against a direct per-bin
recomputation at two test regions (top candidate + a random control
window), exact match. Known limitation, not fixed: a bin at literal 0.0 in
either source track conflates "confirmed closed chromatin" with "no dog
had evidence there" (the same ambiguity already present in the 71-dog
track alone) — a 71:5 blend still applies the full group weight even when
one side is silent, since the two cases aren't distinguishable from the
bigwig values alone.

**76-dog consensus ATAC peaks** (`ATAC_joined76/consensus_peaks_ATAC.bed`):
NOT derived from the joined signal track — recomputed from scratch via
`bedtools multiinter` across all 76 dogs' *individual* MACS3 peak calls
(71 original + 5 new, same per-dog peak files, just combined), then the
same `call_consensus_peaks.py` majority-vote logic as the 71-dog version,
threshold rescaled to **≥ ⌈(76+1)/2⌉ = 39 of 76 dogs (51.3%)**. Result:
9,141 consensus peaks (vs. an unrecorded count for the 71-dog set at the
time — see `05_SHIP/VERSIONS.md` V9-B entry for the downstream effect on
candidate scoring, which turned out to be negligible).

## RRBS methylation (same 71 dogs, PBMC, PRJNA1049514)

Script: `scripts/build_methylation_track.py`. Same two-layer weighting
philosophy as ATAC, adapted to bisulfite data:

**Per-dog global weight**: mapping efficiency and bisulfite conversion rate
(from Bismark), each min-max normalized across the 71 dogs, averaged,
rescaled onto `[0.2, 1.0]` — identical formula shape to the ATAC weight,
different input metrics.

**Per-bin local confidence**: `depth_i(bin) / (depth_i(bin) + λ_i)`, where
λ_i is that dog's own mean CpG depth genome-wide (RRBS is reduced
representation, concentrated near MspI cut sites, so coverage is extremely
uneven — this is what lets an uncovered bin be treated as "not measured"
rather than "0% methylated").

**RRBS_weighted_mean.bw**: same weighted-average formula as ATAC's, but
over %CpG methylation instead of CPM, expressed 0–100.

**RRBS_variability.bw**: standard deviation (not IQR, unlike ATAC — computed
via running sums in a single per-dog streaming pass rather than holding all
71 dogs' arrays in memory at once) of %methylation across dogs with
confidence ≥ 0.2 at that bin.

**RRBS_cpg_coverage_frequency.bw**: per-bin count of dogs with local
confidence ≥ 0.2 at that bin.

## Hi-C (1 dog, "Mischka", Wang et al. 2021, PRJNA587469)

Script: `scripts/build_hic_tracks.py`. Not a population track — a single
German Shepherd's 3 Hi-C libraries, merged, used as structural context only
(explicitly down-weighted in scoring, see `05_SHIP/README.md`).

Pipeline: `pairtools merge` (3 libraries → 1 valid-pairs file) →
`cooler cload pairs` (1 kb fixed-width bins — no restriction-fragment
digest, since the Dovetail kit's enzyme was never disclosed in the source
paper) → `cooler zoomify` (4DN multi-resolution ladder, ICE-balanced at
every zoom level) → `cooltools insulation` at **25 kb resolution**, computed
independently at **window sizes 100 kb, 250 kb, and 500 kb**, boundary
threshold method **"Li"** (cooltools' built-in Li minimum-cross-entropy
thresholding).

**HiC_tad_insulation.tsv / HiC_tad_insulation_{100000,250000,500000}.bw**:
raw insulation score per window size.

**TAD boundary calls** (used as the hard-veto input): a bin is called a TAD
boundary if **any** of the three window sizes flags it as a boundary
(`is_boundary_<window>` from cooltools, OR'd across the three).

## Bin sizes

- ATAC / RRBS Mother Tracks: genome-wide bins as configured in
  `config.yaml` / `config_rrbs.yaml` (`params.bin_size`).
- Hi-C contact matrix base resolution: 1 kb; insulation computed at 25 kb.

## Why this matters for interpretation

Every "weighted mean," "consensus," and "boundary" value in this dataset is
a population-aware, confidence-weighted statistic across up to 71 dogs —
never a single animal's raw signal, and never a naive average that would
silently treat "not measured in this dog" as "closed/unmethylated in this
dog." Anyone re-deriving vetoes or scores from these tracks should use the
formulas above, not assume a generic ATAC/RRBS-seq pipeline.
