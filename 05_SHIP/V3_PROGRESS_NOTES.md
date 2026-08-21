# V3 (ultraconserved elements) — overnight progress notes, 2026-08-21

## Done
- HAL download: 865GB, 100% complete.
- Confirmed dog genome name in HAL: `Canis_lupus_familiaris`. Confirmed it is
  **CanFam3.1**-referenced inside this HAL (chr1-38+X, chrUn_JH373xxx/AAEX0302xxx
  scaffolds), NOT ROS_Cfam_1.0. This required an extra liftover step.
- Reverse liftOver chain found and confirmed working: UCSC
  `GCF_014441545.1ToCanFam3.over.chain.gz` (ROS_Cfam_1.0 -> CanFam3.1).
  All 43 V2 candidates lifted 43/43, zero unmapped. Result:
  `/mnt/d/Jin2024_work/zoonomia_hal/liftover/candidates_canFam3.bed`
  (ordering: line N corresponds to `cand_N` in the MAF filenames below —
  verify this before trusting it downstream, it was assumed not re-checked
  after the fact).
- **Pivoted away from genome-wide hal2maf** after a single 10Mb test region
  produced 11.8GB and hadn't finished after 90+ minutes (whole genome would
  have been ~2.7TB and days). Realized we only need conservation AT the
  candidate loci, not genome-wide — hal2maf on a single candidate-sized
  (~58kb) region took only 53 seconds. This is the single most important
  methodological decision of the night; document it clearly in any writeup.
- hal2maf batch complete: 43 candidate MAFs + 30 random 100kb neutral-model
  regions (genome-wide, chromosome-length-weighted random sampling, seed=42),
  all dog(CanFam3.1)-referenced, `--onlyOrthologs --noAncestors`, in
  `/mnt/d/Jin2024_work/zoonomia_hal/maf_dog/` (14GB total, 73 files).
- **Deliberate simplification vs. the real Zoonomia production pipeline**:
  the actual `modelGeneration_mammals` script filters neutral-model input to
  ancestral-repeat regions specifically (RepeatMasker + halLiftover to
  ancestral-node coordinates + cross-branch synteny validation) before
  phyloFit, to avoid contaminating the neutral model with purifying-selected
  sequence. We are using the 30 random 100kb regions **unfiltered** —
  whole-region MAFs, not restricted to repeat-masked columns. This is a
  reasonable but real simplification (many published methods use "random
  intergenic" or "4-fold-degenerate-site" proxies without the full
  ancestral-repeat-synteny apparatus) — **must be stated explicitly** wherever
  these V3 results get written up or used. Not a silent shortcut.
- score_ship_candidates.py: added `--repeat-content-tsv`/`--repeat-content-threshold`
  (RepeatMasker veto), reran as V6: 43 -> 40 survivors (3 excluded, all >50%
  repeat content, LINE/L1-dominated). Candidate #1 unchanged
  (`NC_051811.1:48,020,921-48,077,046`). File: `candidates_scored_v6.tsv`
  (never overwrite v5 or earlier).
- Consensus-ATAC-peak threshold validated against two independent, dog-native
  regulatory proxies (TSS+-2kb from GFF3, and CpG islands computed natively
  on ROS_Cfam_1.0 via Gardiner-Garden/Frommer criteria — NOT cross-species
  liftover, deliberately, since enhancer liftover across species is
  unreliable): 10.9x enrichment at TSS, 15.1x at CpG islands. Real support
  for keeping the veto as-is. Written up: `05_SHIP/consensus_peak_validation.md`.
- Dog10K SV VCF downloaded (Manta-SV, 1,879 dogs, UU_Cfam_GSD_1.0-referenced).
- UU_Cfam_GSD_1.0 genome downloaded (GCF_011100685.1, 2.5GB) for the
  Dog10K SV liftover-by-realignment step (chain-file-free approach — bwa mem
  realignment of just the survivor candidates, same technique used earlier
  for the original chr12 candidate liftover — NOT a rigorous genome-wide
  liftover, say so wherever this gets used).

## WSL crash and recovery (2026-08-21, ~00:40-01:00)
The first phyloFit attempt concatenated all 30 neutral regions into one file
(7.76GB, 129M lines) and fed that to `phyloFit --subst-mod REV --EM`. This
appears to have driven the WSL2 VM out of memory and crashed the whole VM —
**both** that phyloFit run and the concurrently-running `bwa index` job died
at the same time with no clean error (`bwa index` was even further along,
missing only the final `.sa` file). This was not two independent failures;
one memory spike likely took down the shared WSL2 instance.

Recovery: `bwa index` was rebuilt from scratch and succeeded fully
(`.bwt`/`.sa`/`.pac`/`.ann`/`.amb` all present). phyloFit was retried on a
**single** ~100kb neutral region instead of 30 concatenated — this is
actually closer to the real Zoonomia `modelGeneration_mammals` script's own
scale (it fits phyloFit on one `$maf_100kb` file, not many concatenated) —
and completed cleanly in 15 minutes real time, memory never exceeded ~1.4GB.
Model: `neutral_model_1region.mod` (kept in `/mnt/d/Jin2024_work/zoonomia_hal/`,
not committed to the repo - regenerable from the pipeline above).

**Honest caveat carried forward**: the neutral model is fit on ONE 100kb
region (not the 30 originally planned, and not ancestral-repeat-filtered
like the real pipeline) — this is a weaker, noisier neutral model than a
production Zoonomia release would use. Treat V3 results as a genuine but
provisional pilot, not a final, publication-grade conservation track.

## phyloP results

Ran `phyloP -i MAF --method LRT --mode CONACC --wig-scores` on all 43
candidate MAFs against `neutral_model_1region.mod` (43/43 succeeded, only
benign "pruned leaves with no alignment match" warnings on 3 candidates).

**Coordinate mapping verified, not just assumed**: `cand_N` (from MAF
filenames) maps to line N of both `candidates_canFam3.bed` and
`candidates_passing_ranked_v5.bed` (liftOver preserves row order for mapped
features, and all 43 mapped — no rows dropped). Cross-checked by matching
the unique `score` field between the two files for all 43: 0 mismatches.

**First pass (per-base mean/max) was misleading on its own** — nearly every
candidate showed a similar max single-base score (~11.0-11.2). Investigated
before trusting it: this is expected behavior, not a bug — a single fully-
conserved base hits the same theoretical LRT ceiling regardless of which
candidate it's in, since that ceiling is a property of the fitted model, not
the surrounding window. A single conserved base is not what "ultraconserved
element" means in the literature (Bejerano et al. and similar define it as a
**sustained stretch**, not an isolated base) — a naive per-base max would
have wrongly flagged almost every candidate as equally "conserved."

**Corrected analysis**: computed a 50bp rolling mean per candidate and
looked for sustained stretches above phyloP 2.0. Real, differentiated
signal, not noise — see `05_SHIP/phylop_ultraconserved_v3_pilot.tsv`,
sorted by max 50bp rolling mean:

- 3 candidates (`NC_051843.1:59,813,343-59,874,197`,
  `NC_051843.1:69,851,587-69,907,041`, `NC_051806.1:17,055,808-17,110,317`)
  have **zero** 50bp windows above 2.0 anywhere — genuinely clean.
- 6 candidates stood out with strong sustained conservation (max 50bp
  rolling mean 6.67-8.54, vs. the next-highest candidate at 6.19 — a real
  gap): `NC_051827.1:34,198,197-34,269,938` (8.54), `NC_051807.1:44,207,725-
  44,265,967` (7.81), `NC_051815.1:74,432,972-74,494,463` (7.31 — already
  had the lowest V6 soft score, 0.047), `NC_051812.1:5,534,832-5,600,865`
  (7.08), `NC_051805.1:32,197,089-32,249,732` (6.92),
  `NC_051815.1:15,804,115-15,863,613` (6.67).
- Rare `-20.000` floor values exist (phyloP's numerical floor for
  degenerate/gap-heavy sites) — mostly 1-3 occurrences per candidate except
  `NC_051843.1:79,168,554-79,226,339` with 110, an outlier worth a manual
  look at alignment quality in that window if it's ever reconsidered.

## Dog10K SV — done, clean result

`bwa mem` realignment of the 40 V6-passing candidates against
UU_Cfam_GSD_1.0 (approximate/best-effort — single-window realignment, not a
rigorous coordinate-by-coordinate liftover). **40/40 realigned confidently**
(MAPQ>=30, 0 unmapped). Checked against the Dog10K population SV set
(Manta-SV, 1,879 dogs): **0/40 candidates overlap any called structural
variant.** `05_SHIP/dog10k_sv_check_v6.tsv`.

## Decision — 2026-08-21, acted on

User reviewed the phyloP findings and said explicitly: "retira esses
candidatos" (remove those candidates). Integrated as a real veto in
`score_ship_candidates.py` (`--ultraconserved-tsv`/`--ultraconserved-threshold`,
default 6.5 on max 50bp rolling-mean phyloP — the threshold that separates
the 6 flagged candidates from the rest of the pilot's distribution), rerun
as **V7**. Input file: `05_SHIP/ultraconserved_stretch_v5candidates.tsv`
(built from the pilot TSV above). Consensus-peak validation did not change
any veto (informational only, confirmed the existing threshold). Dog10K SV
did not change the survivor count either way (0/40 had any hit) — not wired
in as a formal veto input since it excluded nobody, documented here instead.

This closes 8/8 of the Ahmed et al. 2026 checklist criteria for the first
time (miRNA, risk-gene radius, gene-dense neighborhood, lncRNA/smallRNA, TAD
content, repeat content, ultraconserved elements, plus the original
TAD-boundary/ATAC-peak/flanking-risk-gene/mappability set) — with the
ultraconserved-elements criterion resting on the single-region-neutral-model
caveat above, not a production-grade Zoonomia conservation track.

## Standing rules being followed
- Never overwrite an earlier `candidates_scored_v*.tsv` / `*.bed` checkpoint.
- Every methodological simplification vs. the "gold standard" pipeline must
  be stated explicitly in writeups, not silently glossed over.
