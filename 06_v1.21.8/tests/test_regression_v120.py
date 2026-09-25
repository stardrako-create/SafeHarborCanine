#!/usr/bin/env python3
"""
Regression tests for the two bug classes fixed 2026-09-10 (cycle 2):

1. keep_bins_with_evidence(): 0 or 1 dog must never be treated as evidence
   of stability - see bw_utils.py's docstring. Both build_mother_track_v2.py
   (ATAC) and build_methylation_track.py (RRBS) call this same function for
   their variability.bw gating, so testing it once covers both.

2. track_value() / build_matched_window_background(): a candidate's own
   value and its background must be read with the same stat_type and
   implicit_zero handling, or the comparison is statistically meaningless
   (the 2026-09-10 bug that made 254/461 candidates spuriously fail the
   accessibility floor). Tested against a small synthetic bigwig with a
   known, hand-computable pattern - if either function silently hardcoded
   a different stat_type instead of respecting the parameter, these
   assertions would catch it.

No pytest dependency - environment.yml has no test framework installed,
and stdlib unittest is enough for both. Run with:
    python3 -m unittest scripts/tests/test_regression_v120.py -v
"""
import os
import sys
import tempfile
import unittest

import numpy as np
import pyBigWig

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from bw_utils import (track_value, build_matched_window_background, keep_bins_with_evidence,
                       evidence_summary)
from score_ship_candidates_v2 import distance_to_nearest_5prime_end


class TestKeepBinsWithEvidence(unittest.TestCase):
    def test_zero_and_one_dog_are_never_evidence(self):
        evidence_count = [0, 1, 2, 3, 1, 0, 5]
        expected = [False, False, True, True, False, False, True]
        got = keep_bins_with_evidence(evidence_count, min_evidence=2)
        np.testing.assert_array_equal(got, expected)

    def test_min_evidence_is_configurable(self):
        evidence_count = [1, 2, 3]
        np.testing.assert_array_equal(
            keep_bins_with_evidence(evidence_count, min_evidence=1),
            [True, True, True],
        )
        np.testing.assert_array_equal(
            keep_bins_with_evidence(evidence_count, min_evidence=3),
            [False, False, True],
        )


class TestMatchedStatistic(unittest.TestCase):
    """A tiny synthetic bigwig, one chromosome, three known-value blocks and
    a trailing gap - big enough (>100_000bp) to pass
    build_matched_window_background's own chrom-size filter."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="shc_regtest_")
        cls.bw_path = os.path.join(cls.tmpdir, "synthetic.bw")
        bw = pyBigWig.open(cls.bw_path, "w")
        bw.addHeader([("chr1", 200_000)])
        # [0, 50000) = 1.0 ; [50000, 100000) = 100.0 ; [100000, 150000) = 5.0
        # [150000, 200000) left as a true gap - no stored data.
        bw.addEntries(
            ["chr1", "chr1", "chr1"],
            [0, 50_000, 100_000],
            ends=[50_000, 100_000, 150_000],
            values=[1.0, 100.0, 5.0],
        )
        bw.close()

    @classmethod
    def tearDownClass(cls):
        try:
            os.remove(cls.bw_path)
            os.rmdir(cls.tmpdir)
        except OSError:
            pass

    def _open(self):
        return pyBigWig.open(self.bw_path)

    def test_stat_type_mean_is_respected(self):
        bw = self._open()
        # equal-width blocks of 1.0 and 100.0 -> bp-weighted mean 50.5
        val = track_value(bw, "chr1", 0, 100_000, stat_type="mean")
        bw.close()
        self.assertAlmostEqual(val, 50.5, places=1)

    def test_stat_type_max_is_respected(self):
        bw = self._open()
        val = track_value(bw, "chr1", 0, 100_000, stat_type="max")
        bw.close()
        self.assertAlmostEqual(val, 100.0, places=4)
        self.assertNotAlmostEqual(val, 50.5, places=1)

    def test_implicit_zero_only_applies_when_requested(self):
        bw = self._open()
        as_missing = track_value(bw, "chr1", 150_000, 200_000, stat_type="mean",
                                  implicit_zero=False)
        as_zero = track_value(bw, "chr1", 150_000, 200_000, stat_type="mean",
                               implicit_zero=True)
        bw.close()
        self.assertIsNone(as_missing)
        self.assertEqual(as_zero, 0.0)

    def test_background_uses_the_same_stat_type_as_the_candidate(self):
        """If build_matched_window_background() ever stopped forwarding
        stat_type to track_value() (e.g. hardcoded "mean" again), a max
        background sampled from this bigwig would cluster far below the
        true stored max (100.0) instead of frequently hitting it."""
        bw = self._open()
        bg = build_matched_window_background(
            bw, window_lengths=[10_000], n_samples=200, seed=1,
            stat_type="max", implicit_zero=True,
        )
        bw.close()
        self.assertEqual(bg.size, 200)
        self.assertGreaterEqual(bg.min(), 0.0)
        self.assertLessEqual(bg.max(), 100.0)
        # a 10kb window landing inside the [50000,100000)=100.0 block
        # should be common enough that the background's own max is high -
        # would fail if "max" silently fell back to a mean-like statistic.
        self.assertGreater(bg.max(), 50.0)


class TestGeneFivePrimeEndDistance(unittest.TestCase):
    """Regression for the 2026-09-11 literature-fidelity fix: Ahmed et al.
    2026 criterion 1 is about distance to the 5' end specifically, which
    depends on strand - a '+' gene's 5' end is its lower coordinate (start),
    a '-' gene's 5' end is its higher coordinate (end). Mixing these up
    would silently reintroduce the exact bug just fixed."""

    def test_plus_strand_5prime_is_start(self):
        genes = {"chr1": [(10_000, 20_000, "GENE_PLUS", "+", 10_000)]}
        # point far from start (5') but close to end (3') - should read as far
        d_near_3prime = distance_to_nearest_5prime_end(genes, "chr1", 19_900)
        d_near_5prime = distance_to_nearest_5prime_end(genes, "chr1", 10_100)
        self.assertEqual(d_near_3prime, 9_900)
        self.assertEqual(d_near_5prime, 100)

    def test_minus_strand_5prime_is_end(self):
        genes = {"chr1": [(10_000, 20_000, "GENE_MINUS", "-", 20_000)]}
        d_near_3prime = distance_to_nearest_5prime_end(genes, "chr1", 10_100)
        d_near_5prime = distance_to_nearest_5prime_end(genes, "chr1", 19_900)
        self.assertEqual(d_near_3prime, 9_900)
        self.assertEqual(d_near_5prime, 100)

    def test_no_gene_exclusion(self):
        """The old (superseded) implementation excluded the two flanking
        genes by name - the corrected one does not, since Ahmed's criterion
        doesn't carve out an exception for whichever genes happen to define
        the candidate window."""
        genes = {"chr1": [(10_000, 10_500, "FLANKING_GENE", "+", 10_000)]}
        d = distance_to_nearest_5prime_end(genes, "chr1", 10_050)
        self.assertEqual(d, 50)  # not infinity - this gene is NOT excluded


class TestEvidenceSummary(unittest.TestCase):
    """Regression for the 2026-09-12 fix (found by Codex, adopted here):
    rrbs_mean's window average must exclude bases a coverage track shows
    were never observed - a dense mean track zero-pads them, which is not
    the same as a measured 0% methylation."""

    @classmethod
    def setUpClass(cls):
        cls.tmpdir = tempfile.mkdtemp(prefix="shc_regtest_evsum_")
        cls.mean_path = os.path.join(cls.tmpdir, "mean.bw")
        cls.cov_path = os.path.join(cls.tmpdir, "cov.bw")
        mean_bw = pyBigWig.open(cls.mean_path, "w")
        mean_bw.addHeader([("chr1", 1000)])
        # Dense, zero-padded convention: covered bases read 80%, everything
        # else (900 of 1000 bases) is written as a literal 0.0, not a gap -
        # non-overlapping, sorted entries, as pyBigWig requires.
        mean_bw.addEntries(["chr1", "chr1", "chr1"], [0, 100, 200],
                            ends=[100, 200, 1000], values=[0.0, 80.0, 0.0])
        mean_bw.close()
        cov_bw = pyBigWig.open(cls.cov_path, "w")
        cov_bw.addHeader([("chr1", 1000)])
        cov_bw.addEntries(["chr1", "chr1", "chr1"], [0, 100, 200],
                           ends=[100, 200, 1000], values=[0.0, 3.0, 0.0])  # 3 dogs covered [100,200)
        cov_bw.close()

    @classmethod
    def tearDownClass(cls):
        for p in (cls.mean_path, cls.cov_path):
            try:
                os.remove(p)
            except OSError:
                pass
        try:
            os.rmdir(cls.tmpdir)
        except OSError:
            pass

    def test_uncovered_bases_are_excluded_not_zero(self):
        mean_bw = pyBigWig.open(self.mean_path)
        cov_bw = pyBigWig.open(self.cov_path)
        mean, frac, n = evidence_summary(mean_bw, cov_bw, "chr1", 0, 1000)
        mean_bw.close()
        cov_bw.close()
        # Naive plain-mean over all 1000 bases (900 zero-padded + 100 real
        # 80%s) would read ~8.0 - the whole point of this fix is that it
        # must NOT do that.
        self.assertAlmostEqual(mean, 80.0, places=4)
        self.assertAlmostEqual(frac, 0.1, places=4)
        self.assertEqual(n, 100)

    def test_zero_coverage_window_returns_none(self):
        mean_bw = pyBigWig.open(self.mean_path)
        cov_bw = pyBigWig.open(self.cov_path)
        mean, frac, n = evidence_summary(mean_bw, cov_bw, "chr1", 300, 400)
        mean_bw.close()
        cov_bw.close()
        self.assertIsNone(mean)
        self.assertEqual(n, 0)


if __name__ == "__main__":
    unittest.main()
