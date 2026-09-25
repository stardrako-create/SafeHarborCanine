from pathlib import Path
import ast,json,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'SV_PBGV000010/alignment_chunked/pair_geometry_20260925'
# Test the geometry function separately from data loading and file generation.
t=ast.parse((A/'review_pair_geometry_20260925.py').read_text(encoding='utf-8-sig'));f=next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name=='geometry');exec(compile(ast.Module(body=[f],type_ignores=[]),'geometry','exec'))
a={'coordinate_target':'chr1','coordinate_start0':100,'coordinate_end0':250,'reverse':False,'AS':150};b={'coordinate_target':'chr1','coordinate_start0':300,'coordinate_end0':450,'reverse':True,'AS':145}
assert geometry(a,b)=={'inward':True,'outer_span_bp':350,'score_sum':295};assert geometry(b,a)==geometry(a,b)
assert geometry(a,dict(b,coordinate_target='chr2')) is None
assert not geometry(dict(a,reverse=True),dict(b,reverse=False))['inward'];assert not geometry(a,dict(b,reverse=False))['inward']
d=json.loads((O/'review.json').read_text());s=[r for r in d['results'] if r['previous_spanning_read_selected']];assert len(s)==7
for r in s:
 assert len(r['sensitivity'])==3
 margins=[x['junction_minus_genomic_score'] for x in r['sensitivity']];assert all(x==margins[0] for x in margins)
 for x in r['sensitivity']:
  assert x['junction']['best_pairs']
  for p in x['junction']['best_pairs']:assert p['inward'] and 197<=p['outer_span_bp']<=360
# No genomic pair means undefined margin, not infinite or conclusive support.
assert sum(r['sensitivity'][0]['junction_minus_genomic_score'] is None for r in s)==2
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'geometry_assertions_passed':5,'spanning_name_groups':7,'with_finite_comparison':5,'without_reported_inward_genomic_alternative':2,'margins_stable_across_tested_bounds':True,'report_sha256':hashlib.sha256((O/'review.json').read_bytes()).hexdigest(),'limitations':'Bounds are descriptive; checks do not calibrate insert-size likelihood or establish independent molecule counts'}
(O/'validation.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
