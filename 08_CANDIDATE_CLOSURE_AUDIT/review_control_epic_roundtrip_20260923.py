from pathlib import Path
import csv,gzip,json,collections,bisect,hashlib,datetime
P=Path(__file__).parent.parent;O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922'
proj=json.loads((O/'EpiC_control_projection.json').read_text())['results'];windows={r['id']:r for r in csv.DictReader((O/'transcript_and_boundary_review.tsv').open(),delimiter='\t')};eligible=[r for r in proj if r['status']=='eligible_unique_contiguous_base_mapping']
def project(points,path):
 by=collections.defaultdict(lambda:collections.defaultdict(list));out={k:set() for k in points}
 for k,values in points.items():
  for chrom,pos in values:by[chrom][pos].append(k)
 coords={c:sorted(v) for c,v in by.items()}
 with gzip.open(path,'rt') as f:
  for line in f:
   x=line.split()
   if not x:continue
   if x[0]=='chain':h=x;tp=int(x[5]);qp=int(x[10]);assert x[4]=='+';continue
   n=int(x[0]);arr=coords.get(h[2],[])
   for pos in arr[bisect.bisect_left(arr,tp):bisect.bisect_left(arr,tp+n)]:
    q=qp+pos-tp;q=q if h[9]=='+' else int(h[8])-1-q
    for key in by[h[2]][pos]:out[key].add((h[7],q))
   tp+=n;qp+=n
   if len(x)==3:tp+=int(x[1]);qp+=int(x[2])
 return out
points={}
for r in eligible:
 for i in range(1000):points[(r['id'],i)]={(r['target'],r['start0']+i if r['strand']=='+' else r['end0']-1-i)}
p1=P/'06_v1.21.7/canFam3ToCanFam6.over.chain.gz';p2=P/'06_v1.21.7/canFam6ToGCF_014441545.1.over.chain.gz'
m1=project(points,p1);m2=project(m1,p2);review=[]
for r in eligible:
 w=windows[r['id']];keys=[(r['id'],i) for i in range(1000)]
 unique1=sum(len(m1[k])==1 for k in keys);unique2=sum(len(m2[k])==1 for k in keys);exact=sum(m2[k]=={(w['chrom'],int(w['start0'])+k[1])} for k in keys)
 review.append(dict(id=r['id'],unique_intermediate_bp=unique1,unique_return_bp=unique2,exact_return_bp=exact,pass_roundtrip=unique1==unique2==exact==1000))
states=json.loads((O/'EpiC_control_states.json').read_text());groups=[]
for r in eligible:
 rs=[a for a in states['annotations'] if a['id']==r['id']];positive=[a for a in rs if any(1<=int(s)<=7 and n>0 for s,n in a['state_bp'].items())]
 groups.append(dict(id=r['id'],promoter_enhancer_tissues=[a['tissue'] for a in positive],states_observed=sorted({int(s) for a in rs for s,n in a['state_bp'].items() if n}),roundtrip_pass=next(x['pass_roundtrip'] for x in review if x['id']==r['id'])))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),roundtrip=review,state_review=groups,chain_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [p1,p2]},state_interpretation='States1-7 promoter/enhancer as defined in existing regulatory_recheck_v1215.py;13 quiescent',limitations=['Coordinate consistency via two additional chain files, not experimental or sequence-identity validation','Same upstream assembly alignment resources may share dependencies','Ten unresolved direct mappings remain unresolved and were not intersected'])
(O/'EpiC_roundtrip_and_state_review.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(tested=len(review),roundtrip_pass=sum(r['pass_roundtrip'] for r in review),state_review=groups),indent=2))
