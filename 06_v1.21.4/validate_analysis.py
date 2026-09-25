from pathlib import Path
import sys,unittest,csv,json,hashlib
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.4'
sys.path.insert(0,str(P/'scripts'))
from build_cpg_islands import independent_calls
class Regression(unittest.TestCase):
 def test_overlap_preserves_separate_coordinates(self):
  d={'CGI_GGF':[(0,1000,1000,.51,.61,120)],'CGI_TJ':[(100,700,600,.60,.70,100)]}
  rows=independent_calls('chr1',d)
  self.assertEqual(len(rows),2);self.assertEqual(rows[0][1:4],(0,1000,'CGI_GGF'));self.assertEqual(rows[1][1:4],(100,700,'CGI_TJ'))
 def test_identical_coordinates_keep_each_definition(self):
  d={t:[(0,600,600,.6,.7,100)] for t in ['CGI_GGF','CGI_TJ']}
  self.assertEqual({r[3] for r in independent_calls('chr1',d)},{'CGI_GGF','CGI_TJ'})
 def test_empty_definition(self):
  self.assertEqual(independent_calls('chr1',{'CGI_TJ':[]}),[])
suite=unittest.defaultTestLoader.loadTestsFromTestCase(Regression)
r=unittest.TextTestRunner(verbosity=2).run(suite)
assert r.wasSuccessful()
with (P/'05_SHIP/cpg_islands_ROS_Cfam_1.0.bed').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
for x in rows:
 n,g,o=(200,.5,.6) if x['type']=='CGI_GGF' else (500,.55,.65)
 assert x['type'] in ('CGI_GGF','CGI_TJ') and int(x['end'])-int(x['start'])==int(x['length'])>=n and float(x['gc_frac'])>=g and float(x['obs_exp'])>=o
assert len({(x['#chrom'],x['start'],x['end'],x['type']) for x in rows})==len(rows)==153908
with (O/'all_local_windows.tsv').open() as f:local=list(csv.DictReader(f,delimiter='\t'))
assert len(local)==2091
for x in local:
 assert int(x['end'])-int(x['start'])==int(x['window_bp'])
 for name in ['repeat_overlap_bp','cpg_overlap_bp','external-regulatory-bed_overlap_bp','atac-peaks-bed_overlap_bp']:
  assert 0<=int(x[name])<=int(x['window_bp'])
 assert 0<=float(x['phylop_scored_fraction'])<=1
chains=json.loads((O/'chain_blocks.json').read_text());best=max(chains['LOC111090199'],key=lambda x:x['covered'])
span=max(b[3] for b in best['blocks'])-min(b[2] for b in best['blocks'])
summary={'cpg_regression_tests':3,'cpg_calls_validated':len(rows),'local_windows_validated':len(local),'loc2_best_chain':{'id':best['id'],'aligned_bp':best['covered'],'target':best['target'],'target_span_bp':span,'aligned_blocks':len(best['blocks']),'target_aligned_fraction':best['covered']/span}}
(O/'validation.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
print('Top ATAC window at each scale:')
for gene in sorted({r['genes'] for r in local}):
 for w in [500,1000,2000]:
  top=max([r for r in local if r['genes']==gene and int(r['window_bp'])==w],key=lambda r:float(r['atac_mean']))
  print(json.dumps(top))
