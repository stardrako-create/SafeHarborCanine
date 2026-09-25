import sys, tempfile, unittest
from unittest.mock import patch
from pathlib import Path
import numpy as np
import pyBigWig
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/06_v1.21.1')
sys.path.insert(0,str(P/'code'))
from bw_utils import track_value,evidence_summary,build_matched_window_background
from build_methylation_track import build_tracks

class CoverageTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory()
  self.handles=[]
  for name,vals in [('mean',[80.,0.,0.,0.]),('cov',[5.,0.,3.,0.])]:
   path=str(Path(self.tmp.name)/(name+'.bw'))
   w=pyBigWig.open(path,'w');w.addHeader([('chr1',200000)])
   w.addEntries(['chr1']*4,[0,50000,100000,150000],ends=[50000,100000,150000,200000],values=vals);w.close()
   self.handles.append(pyBigWig.open(path))
  self.mean,self.cov=self.handles
 def tearDown(self):
  for h in self.handles:h.close()
  self.tmp.cleanup()
 def test_missing_zero_excluded_but_observed_zero_retained(self):
  v,frac,n=evidence_summary(self.mean,self.cov,'chr1',0,200000)
  self.assertEqual(v,40.);self.assertEqual(frac,.5);self.assertEqual(n,100000)
 def test_old_dilution_reproduced_and_corrected(self):
  self.assertEqual(track_value(self.mean,'chr1',0,100000),40.)
  self.assertEqual(track_value(self.mean,'chr1',0,100000,coverage_bw=self.cov),80.)
 def test_all_missing_is_unknown(self):
  self.assertIsNone(track_value(self.mean,'chr1',50000,100000,coverage_bw=self.cov))
 def test_invalid_coordinate_is_not_zero(self):
  with self.assertRaises(ValueError):track_value(self.mean,'wrong',0,10,implicit_zero=True)
 def test_background_is_masked_too(self):
  a=build_matched_window_background(self.mean,[10000],n_samples=25,seed=1,coverage_bw=self.cov)
  self.assertEqual(set(a),{0.,80.})
 def test_empty_background_fails_bounded(self):
  with self.assertRaises(ValueError):build_matched_window_background(self.mean,[10000],n_samples=1,coverage_bw=self.cov,min_evidence=100)
 def test_wrong_mask_operation_rejected(self):
  with self.assertRaises(ValueError):track_value(self.mean,'chr1',0,100,stat_type='max',coverage_bw=self.cov)

class BuilderTests(unittest.TestCase):
 def test_builder_persists_gaps_and_true_measured_zero(self):
  with tempfile.TemporaryDirectory() as td:
   mean,cov,var=[str(Path(td)/f'{n}.bw') for n in ['mean','cov','var']]
   # Two observed bins (one truly unmethylated), two absent bins.
   with patch('build_methylation_track.bin_one_dog',return_value=(np.array([10.,0.,10.,0.]),np.array([8.,0.,0.,0.]))):
    build_tracks(['a','b'],[1.,1.],{'a':{'mean_cpg_depth':1},'b':{'mean_cpg_depth':1}},td,[('chr1',200)],50,mean,cov,var,.2)
   with pyBigWig.open(mean) as m,pyBigWig.open(cov) as c:
    self.assertIsNone(m.stats('chr1',50,100,exact=True)[0])
    self.assertEqual(m.stats('chr1',100,150,exact=True)[0],0.)
    self.assertEqual(track_value(m,'chr1',0,200,coverage_bw=c),40.)

if __name__=='__main__':unittest.main(verbosity=2)
