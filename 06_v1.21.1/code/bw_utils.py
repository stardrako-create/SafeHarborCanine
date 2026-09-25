#!/usr/bin/env python3
"""
Shared helpers for reading bigwig region/window statistics consistently
across the v1.20.0 scripts - written once instead of reimplemented per call
site, which is what caused two real bugs:

1. 2026-09-08: score_ship_candidates.py's peak_frequency background was
   built from only the 2.15% of the genome bedtools multiinter actually
   stores (non-zero-vote regions), silently treating both the background
   and any candidate's own "no data here" case as if they didn't exist,
   rather than as the true zero they are. track_value()'s implicit_zero
   parameter is the fix, applied consistently wherever a sparse
   vote-count/union track is read.

2. 2026-09-10: four of score_ship_candidates_v2.py's five score components
   compared a window-averaged (or window-maxed) candidate value against a
   background built from individual raw bins - a statistical mismatch
   (window aggregation mechanically compresses variance relative to single
   bins, and a window MAX is a biased-high order statistic versus
   single-bin draws regardless). track_value() being the one function both
   a candidate's own value and, via build_matched_window_background(), its
   background sample are computed through means the two sides can no
   longer silently diverge in stat_type - see that function's docstring
   for the 4-condition test that found the original instance of this bug.

keep_bins_with_evidence() is unrelated to bigwig reading: see its own
docstring for the third bug this file fixes (0-or-1-dog "stability").
"""
import random

import numpy as np


def evidence_summary(bw, coverage_bw, chrom, start, end, min_evidence=1):
    """Mean of measured base-pair values, coverage fraction, and contributing bases.
    This is NOT a per-CpG or read-depth-weighted estimator. No coverage is unknown.
    Invalid coordinates/read failures propagate rather than becoming zeros.
    """
    if bw.chroms(chrom) != coverage_bw.chroms(chrom):
        raise ValueError("Mismatched track chromosome lengths")
    v = bw.values(chrom, start, end, numpy=True)
    cov = coverage_bw.values(chrom, start, end, numpy=True)
    mask = np.isfinite(v) & np.isfinite(cov) & (cov >= min_evidence)
    n = int(mask.sum())
    return (float(np.mean(v[mask], dtype=np.float64)) if n else None,
            n / (end-start), n)


def track_value(bw, chrom, start, end, stat_type="mean", implicit_zero=False,
                exact=True, coverage_bw=None, min_evidence=1):
    if chrom not in bw.chroms() or not 0 <= start < end <= bw.chroms(chrom):
        raise ValueError("Invalid genomic interval")
    if coverage_bw is not None:
        if stat_type != "mean" or implicit_zero:
            raise ValueError("Coverage masking currently supports measured means only")
        return evidence_summary(bw, coverage_bw, chrom, start, end, min_evidence)[0]
    val = bw.stats(chrom, start, end, type=stat_type, exact=exact)[0]
    if val is None:
        return 0.0 if implicit_zero else None
    if not np.isfinite(val):
        raise ValueError("Non-finite track statistic")
    return val


def build_matched_window_background(bw, window_lengths, n_samples=3000, seed=42,
                                     stat_type="mean", implicit_zero=False, exact=True,
                                     coverage_bw=None, min_evidence=1):
    """Random genome-wide samples, each a window statistic (stat_type) at a
    length drawn from window_lengths (the real candidate-length
    distribution, ~50-75kb), so the background is the SAME statistic as
    what candidates are compared against - via track_value(), never a
    separately-reimplemented read. Fixes a real bug found 2026-09-10: the
    original accessibility floor compared a window-averaged candidate value
    against a background built from individual raw 25bp bins - an
    apples-to-oranges statistical mismatch (window-averaging mechanically
    compresses variance relative to individual bins, so nearly every window
    mean looks unremarkably close to the middle against that background
    regardless of the window's true position). Isolated with a 4-condition
    test (current production vs. narrowing-only vs. background-matching-
    only vs. both): background-matching alone dropped the failure count
    from 254/461 to 9/461; narrowing alone made it slightly worse (269/461)
    - the bug was the mismatched statistic, not the window width.

    Samples a genome position uniformly (length-weighted chromosome choice),
    not a chromosome uniformly - the earlier version picked a chromosome
    with equal probability regardless of length, over-representing small
    contigs relative to their actual share of the genome."""
    chrom_sizes = dict(bw.chroms())
    chrom_list = [(c, s) for c, s in chrom_sizes.items() if s > 100_000]
    chrom_names = [c for c, _ in chrom_list]
    total_len = sum(s for _, s in chrom_list)
    chrom_weights = [s / total_len for _, s in chrom_list]
    rng = random.Random(seed)
    vals = []
    attempts = 0
    while len(vals) < n_samples:
        attempts += 1
        if attempts > n_samples * 100:
            raise ValueError("Insufficient observable background windows")
        length = rng.choice(window_lengths)
        chrom = rng.choices(chrom_names, weights=chrom_weights, k=1)[0]
        size = chrom_sizes[chrom]
        if size <= length + 2000:
            continue
        pos = rng.randint(1000, size - length - 1000)
        v = track_value(bw, chrom, pos, pos + length, stat_type=stat_type,
                         implicit_zero=implicit_zero, exact=exact,
                         coverage_bw=coverage_bw, min_evidence=min_evidence)
        if v is not None:
            vals.append(v)
    arr = np.array(vals, dtype=np.float64)
    arr.sort()
    return arr


def build_narrow_window_background(bw, radius=10_000, n_samples=3000, seed=42, exact=True):
    """Same matched-statistic idea, for the narrow (likely-insertion-point)
    accessibility check: fixed-radius window centered at random genome
    positions (length-weighted chromosome choice, see
    build_matched_window_background()), radius matching the candidate-side
    narrow-window read."""
    chrom_sizes = dict(bw.chroms())
    chrom_list = [(c, s) for c, s in chrom_sizes.items() if s > 2 * radius + 2000]
    chrom_names = [c for c, _ in chrom_list]
    total_len = sum(s for _, s in chrom_list)
    chrom_weights = [s / total_len for _, s in chrom_list]
    rng = random.Random(seed)
    vals = []
    attempts = 0
    while len(vals) < n_samples:
        attempts += 1
        if attempts > n_samples * 100:
            raise ValueError("Insufficient observable background windows")
        chrom = rng.choices(chrom_names, weights=chrom_weights, k=1)[0]
        size = chrom_sizes[chrom]
        pos = rng.randint(radius + 500, size - radius - 500)
        v = track_value(bw, chrom, pos - radius, pos + radius, stat_type="mean", exact=exact)
        if v is not None:
            vals.append(v)
    arr = np.array(vals, dtype=np.float64)
    arr.sort()
    return arr


def keep_bins_with_evidence(evidence_count, min_evidence=2):
    """Boolean mask: True where at least min_evidence dogs passed the
    per-bin confidence gate. Used to decide which bins get a real entry in
    a variability track and which are left as a genuine bigwig gap - fixes
    a real bug (2026-09-10) where 0-or-1-dog bins wrote variability 0.0,
    read downstream as *maximally stable* - the opposite of *no/
    insufficient evidence*: with 0 dogs, an IQR/stdev over an empty slice
    is NaN, silently nan_to_num'd to 0.0; with 1 dog, the IQR/stdev of a
    single value really is 0 but is equally uninformative about stability
    across dogs. min_evidence=2 is the minimum needed for a spread
    statistic to mean anything at all."""
    return np.asarray(evidence_count) >= min_evidence
