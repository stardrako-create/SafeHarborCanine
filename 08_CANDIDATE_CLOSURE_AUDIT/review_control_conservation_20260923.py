from pathlib import Path
import json,datetime,subprocess,traceback,hashlib,math
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=C/'conservation_20260923';O.mkdir(exist_ok=False)
Z=Path('/mnt/d/Jin2024_work/zoonomia_hal');B='/home/stardrako/miniforge3/envs/cactus/bin/';model=A.parent/'06_v1.21.9/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save():
 t=O/'status.tmp';t.write_text(json.dumps(status,indent=2));t.replace(O/'status.json')
ids={r['id'] for r in json.loads((C/'EpiC_roundtrip_and_state_review.json').read_text())['roundtrip'] if r['pass_roundtrip']}
windows=[r for r in json.loads((C/'EpiC_control_projection.json').read_text())['results'] if r['id'] in ids];assert len(windows)==9
status={'state':'running','started_utc':now(),'completed':[],'errors':[],'model_sha256':hashlib.sha256(model.read_bytes()).hexdigest(),'limitations':['CanFam3-aligned published-model scores transferred only through strict roundtrip windows','HAL extraction suppresses paralogous duplicates; not a separate paralogy assessment','No missing score is interpreted as zero conservation','Remaining 11 controls require resolved coordinate projection']};save()
for w in windows:
 key=w['id'];maf=O/(key+'.maf');wig=O/(key+'.wig');status['active']=key;save()
 try:
  cmd=[B+'hal2maf','--refGenome','Canis_lupus_familiaris','--refSequence',w['target'],'--start',str(w['start0']),'--length','1000','--noAncestors','--noDupes',str(Z/'241-mammalian-2020v2.hal'),str(maf)]
  (O/(key+'.command.json')).write_text(json.dumps(cmd))
  with (O/(key+'.hal.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True,timeout=600)
  positions=set()
  for line in maf.read_text().splitlines():
   s=line.split()
   if len(s)==7 and s[0]=='s' and s[1]=='Canis_lupus_familiaris.'+w['target']:
    assert s[4]=='+';p=int(s[2]);n=int(s[3]);assert len(s[6].replace('-',''))==n
    span=set(range(p,p+n));assert not positions.intersection(span);positions.update(span)
  assert positions==set(range(w['start0'],w['end0'])),('MAF coverage',len(positions))
  cmd=[B+'phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(maf)]
  with wig.open('w') as out,(O/(key+'.phyloP.log')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
  vals={};p=None;step=1
  for line in wig.read_text().splitlines():
   if line.startswith('fixedStep'):
    f=dict(x.split('=') for x in line.split()[1:]);p=int(f['start'])-1;step=int(f.get('step',1));assert int(f.get('span',1))==1
   elif line.strip() and not line.startswith(('#','track')):
    assert p is not None and p not in vals;v=float(line);assert math.isfinite(v);vals[p]=v;p+=step
  assert set(vals)==positions,('WIG coordinates',len(vals),min(vals),w['start0'])
  v=[vals[i] for i in sorted(vals)]
  status['completed'].append({'id':key,'scored_bp':len(v),'mean':sum(v)/len(v),'max50bp_mean':max(sum(v[i:i+50])/50 for i in range(len(v)-49)),'fraction_ge_2':sum(x>=2 for x in v)/len(v),'maf_sha256':hashlib.sha256(maf.read_bytes()).hexdigest(),'wig_sha256':hashlib.sha256(wig.read_bytes()).hexdigest()})
 except Exception as e:status['errors'].append({'id':key,'error':str(e)});(O/(key+'.error.txt')).write_text(traceback.format_exc())
 status['updated_utc']=now();save()
status['state']='completed' if not status['errors'] else 'completed_with_unresolved_errors';status.pop('active',None);save();print(json.dumps(status))
