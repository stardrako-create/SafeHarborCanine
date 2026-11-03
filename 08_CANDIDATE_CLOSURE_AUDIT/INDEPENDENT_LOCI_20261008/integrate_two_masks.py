from pathlib import Path
import json,csv,collections,hashlib
N=Path(__file__).resolve().parent
ids={'NC_051814.1_43352999_43353999','NC_051807.1_2091307_2092307'}
rows=json.loads((N/'normalized_complete_20261008/review.json').read_text())['results'];maps={r['id']:r for r in json.loads((N/'basewise_projection.json').read_text())['results']};reb=json.loads((N/'SNP_reference_rebase_counts_20261008.json').read_text())['records']
out=[]
for w in sorted(ids):
 start=int(w.rsplit('_',2)[1]);mask=set();reasons=[]
 for m in maps[w]['mismatches']:mask.add(m['ROS_offset0'])
 mask.update(maps[w]['ROS_inserted_offsets'])
 for r in rows:
  if r['window']==w and r['source_allele']['filters']==['PASS']:
   lo=r['start0']-start;hi=lo+len(r['ref']);mask.update(range(max(0,lo),min(1000,hi)))
 for r in reb:
  if r['window']==w and r['filters']==['PASS']:mask.add(r['ROS_offset0'])
 counts=collections.defaultdict(lambda:dict(total=0,unique=0,unique_outside_REF_mask=0))
 with (N/'exact_specificity_next9/all_short_segments_exact_uniqueness.tsv').open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['window_id']!=w:continue
   k=int(r['length']);c=counts[k];c['total']+=1;unique=r['unique_exact']=='True';c['unique']+=unique
   lo=int(r['start0'])-start;hi=int(r['end0'])-start;c['unique_outside_REF_mask']+=unique and not any(i in mask for i in range(lo,hi))
 out.append(dict(window=w,masked_bases=len(mask),mask_offsets=sorted(mask),segments=dict(counts),approved=False))
(N/'two_loci_REF_mask_specificity.json').write_text(json.dumps(dict(results=out,limitations=['Mask is observed PASS REF spans and assembly differences, not all indel effects.','No mismatch/bulge/PAM off-target assessment.','Conservation and regional ATAC still pending.']),indent=2));print(json.dumps([{k:v for k,v in r.items() if k!='mask_offsets'} for r in out],indent=2))
