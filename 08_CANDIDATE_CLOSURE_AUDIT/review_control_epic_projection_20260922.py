from pathlib import Path
import csv,gzip,json,collections,hashlib,datetime
P=Path(__file__).parent.parent;O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922';src=Path('/mnt/d/Jin2024_work/zoonomia_hal/liftover/GCF_014441545.1ToCanFam3.over.chain.gz')
windows=list(csv.DictReader((O/'transcript_and_boundary_review.tsv').open(),delimiter='\t'));by=collections.defaultdict(list);maps={}
for w in windows:
 by[w['chrom']].append(w);maps[w['id']]=[set() for _ in range(1000)]
with gzip.open(src,'rt') as f:
 for line in f:
  x=line.split()
  if not x:continue
  if x[0]=='chain':
   h=x;tp=int(x[5]);qp=int(x[10]);assert x[4]=='+';continue
  n=int(x[0])
  for w in by.get(h[2],[]):
   s=int(w['start0']);e=int(w['end0']);lo=max(tp,s);hi=min(tp+n,e)
   for p in range(lo,hi):
    q=qp+p-tp;q=q if h[9]=='+' else int(h[8])-1-q
    maps[w['id']][p-s].add((h[7],q,h[9]))
  tp+=n;qp+=n
  if len(x)==3:tp+=int(x[1]);qp+=int(x[2])
projection=[];eligible={}
for w in windows:
 m=maps[w['id']];unmapped=sum(not x for x in m);ambiguous=sum(len(x)>1 for x in m);r=dict(id=w['id'],unmapped_bp=unmapped,multiple_target_bp=ambiguous)
 if not unmapped and not ambiguous:
  v=[next(iter(x)) for x in m];step=1 if v[0][2]=='+' else -1
  contiguous=all(x[0]==v[0][0] and x[2]==v[0][2] and x[1]==v[0][1]+i*step for i,x in enumerate(v))
  r.update(contiguous=contiguous,target=v[0][0],strand=v[0][2],start0=min(x[1] for x in v),end0=max(x[1] for x in v)+1)
  if contiguous:eligible[w['id']]=r
 r['status']='eligible_unique_contiguous_base_mapping' if w['id'] in eligible else 'unresolved_projection';projection.append(r)
(O/'EpiC_control_projection.json').write_text(json.dumps(dict(chain_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),results=projection,limitations=['Uniqueness within supplied chain set only; not independent sequence or reciprocal validation','Only fully mapped single-target contiguous 1kb windows intersected; unresolved cases not cleared']),indent=2))
annotations=[]
for p in sorted((P/'06_v1.21.5/external_sources').glob('*_13_dense.bed.gz')):
 tissue=p.name.split('_')[0];hits=collections.defaultdict(list)
 with gzip.open(p,'rt') as f:
  for line in f:
   if line.startswith(('#','track','browser')):continue
   x=line.split()
   if len(x)<4:continue
   s=int(x[1]);e=int(x[2]);state=int(x[3])
   for key,r in eligible.items():
    if x[0]==r['target']:
     lo=max(s,r['start0']);hi=min(e,r['end0'])
     if lo<hi:hits[key].append((lo,hi,state))
 for key,r in eligible.items():
  intervals=hits[key];covered=set();counts=collections.Counter()
  for s,e,state in intervals:
   covered.update(range(s,e));counts[state]+=e-s
  annotations.append(dict(id=key,tissue=tissue,covered_bp=len(covered),sum_state_bp=sum(counts.values()),state_bp=dict(counts),nonquiescent_bp=sum(n for s,n in counts.items() if s!=13)))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),windows=len(windows),strictly_projected=len(eligible),unresolved_ids=[r['id'] for r in projection if r['status']=='unresolved_projection'],window_tissue_records=len(annotations),nonquiescent_ids=sorted({r['id'] for r in annotations if r['nonquiescent_bp']}),incomplete_coverage_ids=sorted({r['id'] for r in annotations if r['covered_bp']!=1000 or r['sum_state_bp']!=1000}),annotations=annotations,limitations=['States reported as source labels; state13 treated as quiescent consistently with existing pipeline','Nonquiescent includes diverse states and is not automatically promoter/enhancer','No CAR-T-specific regulatory clearance; independent mapping validation remains pending'])
(O/'EpiC_control_states.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='annotations'},indent=2))
