from pathlib import Path
import json,hashlib,datetime
A=Path(__file__).resolve().parent;O=A/'CONTROL_LOCUS_REVIEW_20260922/repeatmasker_context_20260923_v2'
d=json.loads((O/'review.json').read_text());c=json.loads((O/'contexts.json').read_text())
def fasta(p):
 out={};name=None
 for line in p.read_text().splitlines():
  if line.startswith('>'):name=line[1:].split()[0];out[name]=''
  else:out[name]+=line.strip()
 return out
raw=fasta(O/'controls_context.fa');masked=fasta(O/'controls_context.fa.masked');assert set(raw)==set(masked)==set(c)
checks=[]
for r in d['results']:
 k=r['id'];assert len(raw[k])==len(masked[k]);assert hashlib.sha256(raw[k].encode()).hexdigest()==c[k]['sha256']
 s=c[k]['window_start0']-c[k]['context_start0'];e=c[k]['window_end0']-c[k]['context_start0']
 n=sum(b=='N' and a!='N' for a,b in zip(raw[k][s:e],masked[k][s:e]));assert n==r['masked_bp_union'],(k,n,r['masked_bp_union']);checks.append({'id':k,'masked_bp':n})
lib=Path('/home/stardrako/miniforge3/envs/cactus/share/RepeatMasker/Libraries/CONS-Dfam_4.0/dog')
hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in lib.iterdir() if f.is_file() and f.name in ['speciesMeta.pm','refinelib','cutlib','longlib','shortlib','mirlib','mirslib','retrolib','sinecutlib','shortcutlib']}
(O/'validation.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':'20 input sequences hash-verified; clipped union annotation matches independently counted masked Ns inside all windows','results':checks,'cached_dog_library_sha256':hashes},indent=2));print(json.dumps(checks))
