from pathlib import Path
import json,csv,collections,hashlib
O=Path(__file__).parent
review=json.loads((O/'exact_specificity_20261006/exact_sequence_review.json').read_text());assert review['status']=='completed'
rows=json.loads((O/'UU_projection/inputs.json').read_text())['windows']
maps={x['id']:x for x in json.loads((O/'basewise_projection.json').read_text())['results']}
alleles=json.loads((O/'allele_projection_20261006T115816Z/review.json').read_text())['alleles']
normalized=json.loads((O/'normalized_exact_blocks_20261006_v2/review.json').read_text())['results']
masks={};out=[]
for w in rows:
 key=f"{w['chrom']}_{w['start']}_{w['end']}";m=maps[key];start=int(w['start']);mask=set(range(1000));covered=set()
 for b in m['exact_blocks']:covered.update(range(b['ROS_offset_start0'],b['ROS_offset_end0']))
 mask-=covered
 # Mask projected original REF spans including symbols and unresolved records.
 for x in alleles:
  if x['window']!=key or x['allele']['filters']!=['PASS']:continue
  a=x['allele'];s=a['pos1']-1;e=s+len(a['ref'])
  for b in m['exact_blocks']:
   lo=max(s,b['UU_start0']);hi=min(e,b['UU_end0'])
   for pos in range(lo,hi):mask.add(b['ROS_offset_start0']+pos-b['UU_start0'] if m['strand']=='+' else b['ROS_offset_end0']-1-(pos-b['UU_start0']))
 for x in normalized:
  if x['window']==key and x['source_allele']['filters']==['PASS']:
   mask.update(range(max(0,x['start0']-start),min(1000,x['start0']-start+len(x['ref']))))
 masks[key]=mask
counts=collections.defaultdict(lambda:collections.Counter())
with (O/'exact_specificity_20261006/all_short_segments_exact_uniqueness.tsv').open() as f:
 for x in csv.DictReader(f,delimiter='\t'):
  key=x['window_id'];k=int(x['length']);start=int(key.rsplit('_',2)[1]);s=int(x['start0'])-start;e=int(x['end0'])-start
  counts[key][f'total_{k}']+=1
  if x['unique_exact']=='True':
   counts[key][f'unique_{k}']+=1
   if not any(i in masks[key] for i in range(s,e)):counts[key][f'unique_REF_mask_clear_{k}']+=1
for key,m in masks.items():out.append(dict(id=key,masked_bases=len(m),masked_offsets=sorted(m),counts=dict(counts[key])))
assert len(out)==64 and sum(sum(v for k,v in x['counts'].items() if k.startswith('total_')) for x in out)==243776
report=dict(results=out,limitations=['REF-span mask only; not full indel effect or haplotype clearance','Assembly differences and unmapped bases masked','Unresolved structural/regulatory and biological criteria remain','Exact uniqueness is not mismatch/bulge off-target clearance; no guide selected'])
(O/'masked_specificity_20261008.json').write_text(json.dumps(report,indent=2))
print(json.dumps([x for x in out if x['id']=='NC_051831.1_39708074_39709074']))