from pathlib import Path
import json,hashlib,subprocess,math,datetime,traceback
N=Path(__file__).resolve().parent;S=N/'remote_HAL_verified_20261009_v2'
assert json.loads((S/'status.json').read_text())['state']=='completed','Extraction must complete successfully first'
O=N/'remote_conservation_verified_20261009_v2';O.mkdir(exist_ok=False)
model=N.parent.parent/'06_v1.21.9/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod'
targets={'chr27':(18845873,18846873,'bg1k_1059'),'chr3':(2199848,2200848,'bg1k_1266'),'chr10':(42483814,42484814,'bg1k_1900')};blocks={k:[] for k in targets};positions={k:{} for k in targets};buf=[]
def consume(lines):
 seq=[l.split() for l in lines if l.startswith(('s ', 's\t'))];refs=[v for v in seq if v[1].startswith('Canis_lupus_familiaris.')]
 if not seq:return
 assert len(refs)==1
 r=refs[0];chrom=r[1].split('.',1)[1];assert chrom in targets and r[4]=='+';lo,hi,_=targets[chrom];s=int(r[2]);dna=r[6].replace('-','').upper();assert len(dna)==int(r[3]);assert lo<=s and s+len(dna)<=hi
 assert seq[0]==r;assert len({len(v[6]) for v in seq})==1
 for i,b in enumerate(dna):
  assert s+i not in positions[chrom];positions[chrom][s+i]=b
 blocks[chrom].append('\n'.join(lines)+'\n\n')
for line in (S/'three_regions.maf').read_text().splitlines():
 if not line.strip():consume(buf);buf=[]
 elif not line.startswith('#'):buf.append(line)
consume(buf)
status={'state':'running','source_sha256':hashlib.sha256((S/'three_regions.maf').read_bytes()).hexdigest(),'model_sha256':hashlib.sha256(model.read_bytes()).hexdigest(),'results':[],'errors':[]}
def save():(O/'status.json').write_text(json.dumps(status,indent=2))
save()
for chrom,(lo,hi,name) in targets.items():
 try:
  assert set(positions[chrom])==set(range(lo,hi)),('Incomplete reference',chrom,len(positions[chrom]))
  identity=json.loads((N/'canFam3_reference_identity_20261009'/ (name.removeprefix('bg1k_')+'.response.json')).read_text())
  assert identity['chrom']==chrom and identity['start']==lo and identity['end']==hi
  assert ''.join(positions[chrom][i] for i in range(lo,hi))==identity['dna'].upper(),'HAL reference differs from independently retrieved canFam3 DNA'
  maf=O/(name+'.maf');maf.write_text('##maf version=1\n\n'+''.join(blocks[chrom]));wig=O/(name+'.wig')
  cmd=['/home/stardrako/miniforge3/envs/cactus/bin/phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(maf)];(O/(name+'.command.json')).write_text(json.dumps(cmd))
  with wig.open('w') as out,(O/(name+'.log')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
  vals={};p=None;step=1
  for line in wig.read_text().splitlines():
   if line.startswith('fixedStep'):
    v=dict(x.split('=') for x in line.split()[1:]);p=int(v['start'])-1;step=int(v.get('step',1));assert int(v.get('span',1))==1
   elif line.strip() and not line.startswith(('#','track')):
    assert p is not None and p not in vals;x=float(line);assert math.isfinite(x);vals[p]=x;p+=step
  assert set(vals)==set(positions[chrom]);values=[vals[i] for i in range(lo,hi)];max50=max(sum(values[i:i+50])/50 for i in range(951))
  status['results'].append(dict(id=name,scored_bp=1000,max50bp_mean=max50,exceeds_legacy_6_5=max50>6.5,reference_sha256=hashlib.sha256(''.join(positions[chrom][i] for i in range(lo,hi)).encode()).hexdigest(),maf_sha256=hashlib.sha256(maf.read_bytes()).hexdigest(),wig_sha256=hashlib.sha256(wig.read_bytes()).hexdigest(),limitation='Reference coverage and identity verified against independent UCSC canFam3 sequence; ROS transfer and other gates assessed separately.'))
 except Exception as e:status['errors'].append(dict(id=name,error=str(e)));(O/(name+'.error.txt')).write_text(traceback.format_exc())
 save()
status['state']='completed' if not status['errors'] else 'completed_with_errors';save();print(json.dumps(status))
