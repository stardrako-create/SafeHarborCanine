from pathlib import Path
import json,urllib.parse,hashlib
O=Path(__file__).parent
P=O.parent.parent
rows=json.loads((O/'UU_projection/inputs.json').read_text())['windows']
gff=P/'01_referencia/ROS_Cfam_1.0/genomic.gff'
genes={};features=[]
with gff.open() as f:
 for line in f:
  if line.startswith('#'):continue
  v=line.rstrip().split('\t')
  if len(v)!=9:continue
  a={k:urllib.parse.unquote(x) for item in v[8].split(';') if '=' in item for k,x in [item.split('=',1)]}
  if v[2] in ['gene','pseudogene']:genes[a['ID']]=a.get('Name',a.get('gene',a['ID']))
  elif v[2] not in ['exon','CDS','region']:features.append((v,a))
results=[]
for w in rows:
 s,e=int(w['start']),int(w['end']);hits=[]
 for v,a in features:
  parents=a.get('Parent','').split(',')
  if v[0]!=w['chrom'] or v[6] not in ['+','-'] or not any(x in genes for x in parents):continue
  pos=int(v[3])-1 if v[6]=='+' else int(v[4])-1
  gap=max(pos-e,s-(pos+1),0)
  if gap<=1000000:hits.append(dict(feature=v[2],id=a.get('ID'),gene=';'.join(genes[x] for x in parents if x in genes),position0=pos,gap_bp=gap,strand=v[6]))
 hits.sort(key=lambda x:x['gap_bp'])
 results.append(dict(id=f"{w['chrom']}_{s}_{e}",nearest=hits[0] if hits else None,features=hits))
assert len(results)==64
out=dict(results=results,source=str(gff),sha256=hashlib.sha256(gff.read_bytes()).hexdigest(),limitations=['Annotated RNA 5prime boundaries, not empirical TSS','No exclusion criterion inferred from RNA class or proximity alone'])
(O/'all_RNA_context.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'windows':len(results),'with_nearby_RNA':sum(bool(x['features']) for x in results)}))