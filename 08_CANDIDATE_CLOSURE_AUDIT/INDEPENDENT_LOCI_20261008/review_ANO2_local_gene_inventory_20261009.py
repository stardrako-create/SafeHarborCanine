from pathlib import Path
import json,gzip,urllib.parse,hashlib,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent;sources=[(P/'01_referencia/ROS_Cfam_1.0/genomic.gff','NC_051831.1',39708074,39709074),(N/'ANO2_UU_annotation_crosscheck_20261009/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gff.gz','NC_049248.1',6741835,6742835)];results=[]
for path,c,lo,hi in sources:
 near=[];opener=gzip.open if path.suffix=='.gz' else open
 with opener(path,'rt') as f:
  for line in f:
   if line.startswith('#'):continue
   x=line.rstrip().split('\t')
   if len(x)!=9 or x[0]!=c or x[2]!='gene':continue
   start=int(x[3])-1;end=int(x[4]);d=lo-end+1 if end<=lo else (start-hi+1 if start>=hi else 0)
   if d>=300000:continue
   a={k:urllib.parse.unquote(v) for item in x[8].split(';') if '=' in item for k,v in [item.split('=',1)]};near.append(dict(gene=a.get('gene'),GeneIDs=[item.split(':',1)[1] for item in a.get('Dbxref','').split(',') if item.startswith('GeneID:')],start0=start,end0=end,strand=x[6],distance_bp=d,attributes=a))
 results.append(dict(source=str(path),window=[c,lo,hi],annotated_genes_within300kb=near,source_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),results=results,limits=['Inventory of annotated gene extents, not all unannotated transcription or functional effects.']);p=N/'ANO2_local_gene_inventory_20261009.json';assert not p.exists();p.write_text(json.dumps(result,indent=2));print(json.dumps([dict(window=r['window'],genes=[dict(gene=g['gene'],GeneIDs=g['GeneIDs']) for g in r['annotated_genes_within300kb']]) for r in results]))
