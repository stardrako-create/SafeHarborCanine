from pathlib import Path
import ast,json,csv,gzip,collections,bisect,datetime,hashlib,subprocess,math,traceback
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922'
# Reuse the exact, previously inspected reciprocal-coordinate calculation; do not execute annotation/output section.
source=(A/'review_split_chain_states_20260923.py').read_text(encoding='utf-8-sig');exec(compile(source.split('\nby=collections.defaultdict(lambda:collections.defaultdict(list))')[0],str(A/'review_split_chain_states_20260923.py'),'exec'))
O=C/'split_conservation_20260923';O.mkdir(exist_ok=False)
prior={r['id'] for r in json.loads((C/'EpiC_roundtrip_and_state_review.json').read_text())['roundtrip'] if r['pass_roundtrip']}
points=[{'id':k[0],'ros_relative_pos0':k[1],'target_chrom':next(iter(v))[0],'target_pos0':next(iter(v))[1]} for k,v in good.items() if k[0] not in prior]
(O/'verified_points.json').write_text(json.dumps(points));groups=collections.defaultdict(list)
for p in points:groups[(p['id'],p['target_chrom'])].append(p)
model=P/'06_v1.21.9/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod';B='/home/stardrako/miniforge3/envs/cactus/bin/';hal='/mnt/d/Jin2024_work/zoonomia_hal/241-mammalian-2020v2.hal'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
status={'state':'running','started_utc':now(),'completed':[],'errors':[],'model_sha256':hashlib.sha256(model.read_bytes()).hexdigest()}
def save():
 t=O/'status.tmp';t.write_text(json.dumps(status,indent=2));t.replace(O/'status.json')
save()
for (key,chrom),pts in sorted(groups.items()):
 stem=key+'_'+chrom;maf=O/(stem+'.maf');wig=O/(stem+'.wig');lo=min(p['target_pos0'] for p in pts);hi=max(p['target_pos0'] for p in pts)+1;status['active']=stem;save()
 try:
  assert hi-lo<100000,'Unexpected extraction span; retained pending'
  cmd=[B+'hal2maf','--refGenome','Canis_lupus_familiaris','--refSequence',chrom,'--start',str(lo),'--length',str(hi-lo),'--noAncestors','--noDupes',hal,str(maf)];(O/(stem+'.command.json')).write_text(json.dumps(cmd))
  with (O/(stem+'.hal.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=f,check=True,timeout=600)
  support={};block=[]
  def emit():
   seq=[l.split() for l in block if l.startswith(('s ','s\t'))];dogs=[s for s in seq if s[1]=='Canis_lupus_familiaris.'+chrom]
   if not dogs:return
   assert len(dogs)==1;dog=dogs[0];assert dog[4]=='+';pos=int(dog[2])
   for s in seq:assert len(s[6])==len(dog[6])
   for i,b in enumerate(dog[6]):
    if b=='-':continue
    assert pos not in support;support[pos]=len({s[1].split('.')[0] for s in seq if not s[1].startswith('Canis_lupus_familiaris.') and s[6][i].upper() in 'ACGT'});pos+=1
  for line in maf.read_text().splitlines()+['']:
   if not line.strip():emit();block=[]
   else:block.append(line)
  assert set(support)==set(range(lo,hi))
  vals={}
  if any(support.values()):
   cmd=[B+'phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(maf)]
   with wig.open('w') as f,(O/(stem+'.phyloP.log')).open('w') as err:subprocess.run(cmd,stdout=f,stderr=err,check=True,timeout=600)
   pos=None;step=1
   for line in wig.read_text().splitlines():
    if line.startswith('fixedStep'):
     fields=dict(x.split('=') for x in line.split()[1:]);pos=int(fields['start'])-1;step=int(fields.get('step',1));assert int(fields.get('span',1))==1
    elif line.strip() and not line.startswith(('#','track')):
     assert pos is not None and pos not in vals;v=float(line);assert math.isfinite(v);vals[pos]=v;pos+=step
   assert set(vals)<=set(support)
  mapped=[]
  for p in pts:
   q=p['target_pos0'];mapped.append(dict(p,non_dog_species=support[q],score=vals.get(q) if support[q]>0 else None,reason='no_comparative_ACGT' if support[q]==0 else 'score_available' if q in vals else 'missing_score'))
  (O/(stem+'.mapped_scores.json')).write_text(json.dumps(mapped,indent=2));scored=[p['score'] for p in mapped if p['score'] is not None]
  status['completed'].append({'id':key,'chrom':chrom,'reciprocal_mapped_bp':len(mapped),'supported_scored_bp':len(scored),'mean_supported_only':sum(scored)/len(scored) if scored else None,'fraction_ge2_supported_only':sum(v>=2 for v in scored)/len(scored) if scored else None,'maf_sha256':hashlib.sha256(maf.read_bytes()).hexdigest(),'mapped_sha256':hashlib.sha256((O/(stem+'.mapped_scores.json')).read_bytes()).hexdigest()})
 except Exception as e:status['errors'].append({'id':key,'chrom':chrom,'error':str(e)});(O/(stem+'.error.txt')).write_text(traceback.format_exc())
 status['updated_utc']=now();save()
status.pop('active',None);status['state']='completed' if not status['errors'] else 'completed_with_errors';status['limitations']=['Only exact reciprocal coordinate points reported, not whole window approval','Unsupported positions have null score even when phyloP emits numeric values','No new conservation threshold introduced; means only over supported scored subset','Shared chain sources and noDupes extraction retain their original limitations'];save();print(json.dumps(status,indent=2))
