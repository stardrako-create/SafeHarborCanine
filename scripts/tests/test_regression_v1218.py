import sys
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

class TestWrittenEvidence(unittest.TestCase):
    def test_zero_one_two_dogs_produce_correct_gaps_and_values(self):
        import tempfile, warnings, pyBigWig
        from build_mother_track_v2 import build_accessibility_and_variability
        with tempfile.TemporaryDirectory() as root:
            root=Path(root)
            for name,values in [('a',[0.,2.,2.]),('b',[0.,0.,6.])]:
                d=root/name/'bigwig';d.mkdir(parents=True)
                with pyBigWig.open(str(d/(name+'.cpm.bw')),'w') as bw:
                    bw.addHeader([('chr1',75),('chr2',25)])
                    bw.addEntries(['chr1']*3,[0,25,50],ends=[25,50,75],values=values)
                    bw.addEntries(['chr2'],[0],ends=[25],values=[0.])
            mean=root/'mean.bw';var=root/'var.bw'
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                build_accessibility_and_variability(['a','b'],[1.,1.],np.array([1.,1.]),np.array([1.,1.]),str(root),[('chr1',75),('chr2',25)],25,str(mean),str(var))
            with pyBigWig.open(str(mean)) as bw:
                v=bw.values('chr1',0,75,numpy=True)[::25]
                np.testing.assert_allclose(v,[np.nan,2./3.,4./3.],equal_nan=True,rtol=1e-6)
                self.assertIsNone(bw.stats('chr2',0,25,exact=True)[0])
            with pyBigWig.open(str(var)) as bw:
                np.testing.assert_allclose(bw.values('chr1',0,75,numpy=True)[::25],[np.nan,np.nan,2.],equal_nan=True)

if __name__=='__main__':unittest.main()
