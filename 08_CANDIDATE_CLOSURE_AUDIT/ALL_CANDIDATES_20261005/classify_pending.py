from pathlib import Path
import json,collections,hashlib
O=Path(__file__).parent
src=O/'allele_projection_20261006T115816Z/review.json'
rows=json.loads(src.read_text())['alleles']
maps={r['id']:r for r in json.loads((O/'basewise_projection.json').read_text())['results']}
projections={r['id']:next(h for h in r['hits'] if h['primary']) for r in json.loads((O/'UU_projection/projection_review.json').read_text())['results']}
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results=[]
for x in rows:
 if x['status']!='needs_boundary_or_reference_reconciliation':continue
 a=x['allele'];m=maps[x['window']];h=projections[x['window']];s=a['pos1']-1;e=s+len(a['ref'])
 if a['alt']=='*':reason='spanning_deletion_symbol_requires_upstream_record'
 elif not a['sequence_allele']:reason='other_nonsequence_allele'
 elif not h['start0']<=s<e<=h['end0']:reason='crosses_alignment_envelope'
 else:
  mism=[v for v in m['mismatches'] if s<=v['UU_position0']<e]
  if len(a['ref'])==len(a['alt'])==1 and len(mism)==1:
   reason='catalogue_ALT_matches_ROS_reference' if a['alt'].upper()==mism[0]['query_in_UU_orientation'] else 'SNP_at_reference_difference_other_ALT'
  elif mism:reason='allele_overlaps_reference_mismatch'
  else:reason='allele_crosses_indel_or_exact_block_boundary'
 results.append(dict(x,reason=reason,source_PASS=a['filters']==['PASS']))
assert len(results)==422
summary={}
for w in maps:
 z=[x for x in results if x['window']==w]
 summary[w]=dict(pending=len(z),PASS_pending=sum(x['source_PASS'] for x in z),reasons=dict(collections.Counter(x['reason'] for x in z)))
report=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),counts=dict(collections.Counter(x['reason'] for x in results)),PASS_counts=dict(collections.Counter(x['reason'] for x in results if x['source_PASS'])),windows=summary,observations=results,limitations=['ALT matching ROS is not absence of variation; UU REF may be a nonreference allele in ROS','Star alleles denote spanning deletions and are not simple SNPs','No failed filters converted to PASS; no unresolved calls discarded'])
(O/'pending_classification_20261006.json').write_text(json.dumps(report,indent=2))
print(json.dumps({k:report[k] for k in ['counts','PASS_counts']}))
print(json.dumps(summary['NC_051831.1_39708074_39709074']))