# Accessibility criteria and a cross-pipeline benchmark proposal

Written 2026-09-11, responding to two design-level points from a review of
the v120c recalibration - not bugs, genuine open questions about what the
pipeline's accessibility/peak criteria are actually supposed to mean.

## 1. What distinguishes "accessible enough" from "a regulatory element to avoid"

The pipeline already draws this distinction on two separate axes, but it
has never been written down explicitly. Making it explicit surfaces real
gaps.

**Axis 1 - discrete peak status** (`veto_atac_peak`, `atac_peak_frequency`
/ `score_low_peak_frequency`): is there a MACS3-called peak here, in any of
the 76 dogs? `veto_atac_peak` hard-excludes direct overlap with the
consensus peak set; `score_low_peak_frequency` continuously rewards
candidates far from where peaks cluster across dogs. This is the proxy for
"this looks like a known/likely regulatory element - don't disrupt it."

**Axis 2 - continuous accessibility level** (`veto_low_atac_accessibility`,
the fold-enrichment `atac_mean`/`atac_mean_narrow`): is chromatin generally
open here, independent of whether a discrete peak was ever called? This is
the proxy for "will DNA inserted here actually be transcribable" - a
functional-insert concern, not a safety concern per se.

NC_051833.1:23,112,498-23,172,320 (the candidate that actually dominates
the ranking at a realistic threshold - see the score comparison below) is a
clean illustration of the two axes working as intended: `atac_peak_frequency
= 0.0` (zero dogs ever called a peak here - not a regulatory element by
this definition) while `atac_mean_percentile = 26.5%` (moderately, not
dramatically, above the genome-wide background - open enough without being
a hotspot).

**The real gaps, stated plainly rather than assumed away:**

- Peak-calling is itself a threshold-dependent, imperfect proxy. `peak_
  frequency = 0` means "MACS3 didn't call a peak here in our specific
  76-dog dataset at our specific settings" - not "this is definitely not a
  regulatory element." A real but weak/context-specific element could sit
  below MACS3's calling threshold in every dog we have and still matter
  functionally. This isn't fixable by re-tuning our own pipeline; it's an
  inherent limit of peak-calling as ground truth, worth stating to
  Ehsan/Vasco rather than treating our own peak set as authoritative.

- The two axes can genuinely pull in opposite directions. The safest
  possible signature by Axis 1 (zero peak activity, ever) doesn't require
  any particular value on Axis 2 - a fully quiet, low-accessibility region
  is just as "no peak" as a moderately-open one. The accessibility floor
  exists for a DIFFERENT reason (insert functionality), and currently we
  require BOTH a peak-safety margin AND a minimum accessibility floor
  simultaneously. That's a defensible design choice, but it is a choice,
  with a real tradeoff - it's not obviously "the" correct criterion, and
  should be presented to Ehsan/Vasco as a decision to weigh in on, not an
  already-settled fact.

## 2. Cross-pipeline ATAC benchmark - a concrete proposal, not yet run

Motivated by Ehsan's own "~30" objection and the resolution-mismatch
finding (our window-averaged candidate scores and his apparent peak-level
numbers are almost certainly not the same statistic). Proposed design,
for discussion before either side runs anything:

1. **Same coordinates, both sides.** Score our current shortlist (the 16
   negative-control candidates, ANO2/NTF3, NC_051833.1, and rank #8) AND
   a handful of Ehsan's own candidates, on both pipelines, at the exact
   same intervals.
2. **Same underlying data where possible.** Ideally both pipelines process
   the identical set of BAM files (the shared 76-dog cohort, or whatever
   subset both sides can agree on) rather than comparing across different
   sample sets - eliminates "different data" as a confound.
3. **Include known reproducible peaks as anchors.** A small set of
   well-established open-chromatin regions in canine PBMC ATAC-seq (e.g.
   housekeeping-gene promoters) scored on both pipelines first - if both
   sides agree on where the KNOWN peaks are, that calibrates the
   comparison; if they disagree even there, that's a more fundamental
   methodological gap worth resolving before comparing candidates at all.
4. **Same resolution, explicitly stated.** Report bin-level/peak-summit
   values AND window-averaged values on both sides, labeled as such - the
   "30 vs ~1" gap this whole thread started from is a resolution mismatch,
   not necessarily a disagreement about the underlying biology; making
   resolution explicit on both sides should resolve or properly localize
   that gap.
5. **Standard QC alongside peak height** (ENCODE ATAC-seq QC guidelines):
   signal-to-noise, FRiP, reproducibility across dogs/replicates - not
   just "does this look like a tall peak."

Not run - needs Ehsan's input on his own pipeline's parameters/outputs and
agreement on which BAMs and reference intervals to use.

## Score comparison referenced above (from candidates_scored_v120c.tsv)

| candidate | genes | final_score | percentile | wide/narrow accessibility pctile | peak_frequency |
|---|---|---|---|---|---|
| NC_051831.1:39677324-39739751 | ANO2/NTF3 | 0.4309 | 49.5 | 73.7 / 67.3 | 1/76 |
| NC_051833.1:23112498-23172320 | LOC111093212/LOC111093135 | 0.6789 | 92.2 | 26.5 / 34.5 | 0/76 |

ANO2/NTF3 clears the strict p55 accessibility bar; NC_051833.1 does not,
but scores substantially higher overall (better ATAC stability, better
methylation, better peak-frequency score) and dominates the ranking once
the accessibility bar is set anywhere looser than ~p35. Its own weak point:
RRBS stability is poor (score 0.18, high inter-dog methylation variance at
that locus) - not a clean win either, just a different tradeoff.
