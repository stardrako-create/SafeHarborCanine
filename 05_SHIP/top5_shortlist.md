# Final shortlist — Top 1 + 4 backups

This closes the last item in the project's own "definição de feito" for the
bioinformatics discovery phase: a principal candidate plus at least one
backup, not just a ranked list. All 5 pass every criterion in
`candidates_scored_v8.tsv` (461 SHIP candidates → 34 survivors → top 5
here), including all 8/8 Ahmed et al. 2026 criteria, the RepeatMasker
repeat-content check, the ultraconserved-elements phyloP check, and the
Dog10K structural-variant check. Ranked by `final_score` (genome-wide
percentile composite — see `VERSIONS.md` for the V8 rewrite).

| Rank | Coordinates (ROS_Cfam_1.0) | Length | Left gene | Right gene | Score |
|---|---|---:|---|---|---:|
| 1 | `NC_051805.1:7,072,137-7,132,579` | 60,442 bp | LOC111090579 | LOC100685067 | 0.7649 |
| 2 | `NC_051811.1:48,020,921-48,077,046` | 56,125 bp | RIT2 | LOC119872716 | 0.7576 |
| 3 | `NC_051835.1:23,436,373-23,508,916` | 72,543 bp | LOC111093569 | LOC119867012 | 0.7548 |
| 4 | `NC_051826.1:19,596,764-19,663,367` | 66,603 bp | LOC100687588 | LOC100687653 | 0.7438 |
| 5 | `NC_051805.1:60,400,908-60,462,010` | 61,102 bp | LOC119869937 | LOC111089986 | 0.7407 |

## Why #1, specifically

`NC_051805.1:7,072,137-7,132,579` ranks highest on the composite score, and
none of its individual components are an outlier driving that alone —
worth checking, since a composite score can hide one dominant axis:

- Repeat content 31.25% (well under the 50% veto threshold, close to
  genome average ~35-36%)
- Ultraconserved-elements signal (max 50bp-rolling phyloP): 4.16, clearly
  below the 6.5 veto threshold and below every one of the 6 candidates
  excluded in V7
- Not one of the 3 candidates independently confirmed by Ehsan Valiollahi's
  separate SHIP filtering (see `ehsan_crossvalidation.md`) — that
  independent-convergence evidence currently belongs to
  `NC_051807.1:10,779,807-10,837,318`,
  `NC_051821.1:4,818,014-4,875,944`, and
  `NC_051843.1:10,578,732-10,643,009`, none of which made the top 5 here.
  Worth keeping in mind: this shortlist is ranked by our own composite
  score alone, not by cross-method agreement — the two aren't the same
  question.

## Why 4 backups, not just #1

Scores 2-5 are close to #1 (0.7407-0.7649, a 3.2% spread) — this is a
shortlist of comparably strong candidates, not one clear winner miles
ahead of the rest. If gRNA design or off-target scoring rules out #1 for a
reason not captured by this pipeline (e.g. no usable PAM site, or an
off-target hit specific to that exact sequence), any of #2-#5 is a
reasonable next attempt without dropping back to the full 34.

## What this shortlist does not cover

No experimental validation, no gRNA design, no off-target scoring yet -
see `VERSIONS.md` / `05_SHIP/README.md` "Known gaps" for what's still
ahead before any of these 5 can be called more than a computational
priority. `confidence` in EpiLog terms remains `untested` for all 5.

## Files
- `top5_shortlist.bed` — the 5 coordinates above, ranked, BED format
- `candidates_scored_v8.tsv` — full 461-candidate table these were drawn from
