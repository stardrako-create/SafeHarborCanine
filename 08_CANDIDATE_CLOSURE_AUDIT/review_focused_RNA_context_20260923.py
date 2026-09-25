from pathlib import Path
import json,urllib.parse,hashlib
A=Path(__file__).resolve().parent;O=A/'W01_FOCUSED_REVIEW_20260923';d=json.loads((O/'focused_checks.json').read_text());windows={k:(r['chrom'],int(r['start0']),int(r['end0'])) for k,r in d['controls'].items()};windows['w01']=('NC_051811.1',17255706,17256706);out={k:{'nearest_all_RNA_5prime':None,'overlap_transcripts':[]} for k in windows}
gff=A.parent/'01_referencia/ROS_Cfam_1.0/genomic.gff'
with gff.open() as f:
 for line in f:
  if line.startswith('#'):continue
  x=line.rstrip().split('\t')
  if len(x)!=9 or not ('RNA' in x[2] or x[2]=='transcript') or x[6] not in ('+','-'):continue
  s=int(x[3])-1;e=int(x[4]);t=s if x[6]=='+' else e-1;attrs=dict(v.split('=',1) for v in x[8].split(';') if '=' in v)
  for k,(c,a,b) in windows.items():
   if x[0]!=c:continue
   gap=a-t if t<a else t-b if t>=b else 0
   r={'feature':x[2],'id':attrs.get('ID'),'gene':urllib.parse.unquote(attrs.get('gene','')),'strand':x[6],'annotated_5prime0':t,'gap_bp':gap}
   if out[k]['nearest_all_RNA_5prime'] is None or gap<out[k]['nearest_all_RNA_5prime']['gap_bp']:out[k]['nearest_all_RNA_5prime']=r
   if s<b and a<e:out[k]['overlap_transcripts'].append(r)
assert out['w01']['nearest_all_RNA_5prime']['gap_bp']==32785
(O/'all_RNA_context.json').write_text(json.dumps({'results':out,'gff_sha256':hashlib.sha256(gff.read_bytes()).hexdigest(),'limit':'Annotated 5prime boundaries, not experimental alternative TSS mapping'},indent=2));print(json.dumps(out,indent=2))
