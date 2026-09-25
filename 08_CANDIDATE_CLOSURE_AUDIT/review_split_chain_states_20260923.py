from pathlib import Path
import ast,json,csv,gzip,collections,bisect,datetime,hashlib
A=Path(__file__).resolve().parent;P=A.parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=C/'gapped_projection_20260923'
tree=ast.parse((A/'review_control_epic_roundtrip_20260923.py').read_text(encoding='utf-8-sig'));func=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='project');exec(compile(ast.Module(body=[func],type_ignores=[]),'project','exec'))
rows={r['id']:r for r in csv.DictReader((C/'transcript_and_boundary_review.tsv').open(),delimiter='\t')};review=json.loads((O/'CanFam3_review.json').read_text());points={}
for r in review['results']:
 for run in r['unique_contiguous_runs']:
  step=1 if run['strand']=='+' else -1
  for i in range(run['ros_relative_start0'],run['ros_relative_end0']):points[(r['id'],i)]={(run['chrom'],run['first_target_pos0']+step*(i-run['ros_relative_start0']))}
p1=P/'06_v1.21.7/canFam3ToCanFam6.over.chain.gz';p2=P/'06_v1.21.7/canFam6ToGCF_014441545.1.over.chain.gz';m1=project(points,p1);m2=project(m1,p2)
good={k:v for k,v in points.items() if len(m1[k])==1 and m2[k]=={(rows[k[0]]['chrom'],int(rows[k[0]]['start0'])+k[1])}}
by=collections.defaultdict(lambda:collections.defaultdict(list))
for k,v in good.items():
 chrom,p=next(iter(v));by[chrom][p].append(k)
coords={c:sorted(v) for c,v in by.items()};annotations=[]
for path in sorted((P/'06_v1.21.5/external_sources').glob('*_13_dense.bed.gz')):
 states={}
 with gzip.open(path,'rt') as f:
  for line in f:
   if line.startswith(('#','track','browser')):continue
   x=line.split()
   if len(x)<4:continue
   arr=coords.get(x[0],[]);lo=bisect.bisect_left(arr,int(x[1]));hi=bisect.bisect_left(arr,int(x[2]))
   for p in arr[lo:hi]:
    for k in by[x[0]][p]:assert k not in states;states[k]=int(x[3])
 for key in rows:
  vals=[v for k,v in states.items() if k[0]==key];annotations.append({'id':key,'tissue':path.name.split('_')[0],'covered_bp':len(vals),'state_bp':dict(collections.Counter(vals))})
summary=[]
for key in rows:
 n=sum(k[0]==key for k in good);t=[a for a in annotations if a['id']==key];summary.append({'id':key,'reciprocal_unique_bp':n,'unresolved_bp':1000-n,'all_tissues_cover_projectable_bases':all(a['covered_bp']==n for a in t),'promoter_enhancer_tissues_in_projectable_bases':[a['tissue'] for a in t if any(1<=s<=7 for s in a['state_bp'])]})
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'summary':summary,'annotations':annotations,'chain_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [p1,p2]},'limitations':['Split reciprocal coordinate consistency only, not sequence identity','Partial regulatory annotation cannot clear unprojectable bases or whole windows','Previously fully contiguous criteria unchanged; no unresolved base treated as quiescent']}
(O/'CanFam3_partial_roundtrip_states.json').write_text(json.dumps(result,indent=2));print(json.dumps(summary,indent=2))
