from pathlib import Path
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');S=P/'scripts'
p=S/'score_ship_candidates_v2.py';s=p.read_text(encoding='utf-8');anchor='def percentile_rank(sorted_background, value):';helper='''def missing_score_components(candidate, weights):
    """A partial score cannot supply evidence for a missing enabled component."""
    return [f"score_component:{name}" for name, weight in weights.items()
            if weight > 0 and candidate.get(f"score_{name}") is None]


''';assert anchor in s
if 'def missing_score_components(' not in s:s=s.replace(anchor,helper+anchor)
old='''        for component, weight in weights.items():
            if weight > 0 and c.get(f"score_{component}") is None:
                missing.append(f"score_component:{component}")''';assert old in s;s=s.replace(old,'        missing.extend(missing_score_components(c, weights))');p.write_text(s,encoding='utf-8')
p=S/'build_mother_track_v2.py';s=p.read_text(encoding='utf-8');old='    return arr.reshape(n_bins, bin_size).mean(axis=1)';new='''    # Padding is an implementation detail, not additional genomic bases.
    lengths = np.minimum(bin_size, size - np.arange(n_bins) * bin_size)
    return (arr.reshape(n_bins, bin_size).sum(axis=1, dtype=np.float64) / lengths).astype(np.float32)''';assert old in s;s=s.replace(old,new);p.write_text(s,encoding='utf-8')
test='''import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from build_mother_track_v2 import read_binned
from score_ship_candidates_v2 import missing_score_components

class FakeBW:
    def __init__(self,mode):self.mode=mode
    def chroms(self):return {} if self.mode=='missing' else {'chr1':5}
    def values(self,c,a,b,numpy=True):
        if self.mode=='corrupt' or (self.mode=='whole_fails' and a==0 and b==5):
            raise RuntimeError('unreadable')
        return np.array([2,2,4,4,8],dtype=np.float32)[a:b]

class TestReadFailures(unittest.TestCase):
    def test_irrecoverable_chunk_is_not_zero(self):
        with self.assertRaisesRegex(RuntimeError,'zero substitution forbidden'):
            read_binned(FakeBW('corrupt'),'chr1',5,2,3,chunk_size=2)
    def test_missing_chromosome_is_not_zero(self):
        with self.assertRaises(RuntimeError):read_binned(FakeBW('missing'),'chr1',5,2,3)
    def test_recovered_chunks_and_partial_bin(self):
        np.testing.assert_array_equal(read_binned(FakeBW('whole_fails'),'chr1',5,2,3,chunk_size=2),[2,4,8])

class TestPartialScore(unittest.TestCase):
    def test_partial_score_still_reports_missing_component(self):
        self.assertEqual(missing_score_components({'final_score':0.9,'score_a':0.9},{'a':1,'b':1}),['score_component:b'])
    def test_zero_value_is_observed(self):
        self.assertEqual(missing_score_components({'score_a':0.0},{'a':1}),[])
    def test_disabled_component_is_not_required(self):
        self.assertEqual(missing_score_components({}, {'a':0}),[])
if __name__=='__main__':unittest.main()
'''
(S/'tests/test_regression_v1218.py').write_text(test,encoding='utf-8');print('Added specific regressions and corrected terminal-bin denominator.')
