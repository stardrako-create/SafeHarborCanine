from pathlib import Path
import json,collections,datetime
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';Q=C/'variant_catalogue_query_20260923';O=C/'star_context_20260923';O.mkdir(exist_ok=False)
review=json.loads((C/'variant_projection_20260923/review.json').read_text());deletions=[]
for f in Q.glob('bg*.json'):
 d=json.loads(f.read_text())
 for n,a in enumerate(d['alleles']):
  ref=a['ref'];alt=a['alt']
  if not a['sequence_allele'] or len(ref)<=len(alt):continue
  left=0
  while left<min(len(ref),len(alt)) and ref[left]==alt[left]:left+=1
  right=0
  while right<min(len(ref)-left,len(alt)-left) and ref[-right-1]==alt[-right-1]:right+=1
  remaining_alt=alt[left:len(alt)-right if right else len(alt)]
  if remaining_alt:continue # complex replacement is not a proven pure deletion
  start=a['pos1']-1+left;end=a['pos1']-1+len(ref)-right
  deletions.append({'source':f.name,'allele_index0':n,'chrom':a['uu_chrom'],'deleted_start0':start,'deleted_end0':end,'allele':a})
result=[]
for p in review['pending_alleles']:
 a=p['allele'];assert a['alt']=='*'
 hits=[d for d in deletions if d['chrom']==a['uu_chrom'] and d['deleted_start0']<=a['pos1']-1<d['deleted_end0']]
 distinct={(h['chrom'],h['deleted_start0'],h['deleted_end0']) for h in hits}
 result.append({'id':p['id'],'position1':a['pos1'],'chrom':a['uu_chrom'],'compatible_local_deletion_records':hits,'distinct_deleted_intervals':len(distinct),'status':'compatible_local_deletion_found_unphased' if hits else 'origin_unresolved_in_downloaded_windows'})
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'star_entries':len(result),'counts':dict(collections.Counter(r['status'] for r in result)),'results':result,'limitations':['Coordinate compatibility only: no sample genotype or phase establishes origin','Source records retain duplicates across catalogues','Only downloaded regional records examined; upstream deletions outside these records can be missed','Do not count star entries as independent deletions or project star as sequence']}
(O/'review.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary['counts']))
