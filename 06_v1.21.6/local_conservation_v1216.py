from pathlib import Path
import csv,subprocess,json,re,numpy as np
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.6';Z=Path('/mnt/d/Jin2024_work/zoonomia_hal')
files={'ANO2':'maf_dog/cand_36.chr27.39320764.maf','LOC119876429':'new3_phylop/loc1_LOC119876429_LOC119872513.maf','NPNT':'new3_phylop/loc3_NPNT_TBCK.maf'}
models={'pilot':Z/'neutral_model_1region.mod','scale_check':O/'independent_region_scale.mod'}
rows=[]
for w in csv.DictReader((P/'06_v1.21.5/local_canfam3_blocks.tsv').open(),delimiter='\t'):
 g=w['genes'].split('/')[0];a=int(w['start']);b=int(w['end']);target=O/(g+'.local.maf');covered=0
 with target.open('w') as out:
  out.write('##maf version=1\n\n');block=[]
  def emit(block):
   global covered
   seq=[l.split() for l in block if l.split() and l.split()[0]=='s'];dog=next((s for s in seq if s[1]=='Canis_lupus_familiaris.'+w['chrom']),None)
   if dog is None:return
   assert dog[4]=='+';start=int(dog[2]);end=start+int(dog[3])
   if end<=a or start>=b:return
   ix=[];pos=start
   for i,c in enumerate(dog[6]):
    if c!='-':
     if a<=pos<b:ix.append(i)
     pos+=1
   if not ix:return
   lo=ix[0];hi=ix[-1]+1;covered+=len(ix);out.write('a score=0\n')
   for s in [dog]+[x for x in seq if x is not dog]:
    text=s[6][lo:hi];size=len(text.replace('-',''))
    if not size:continue
    newstart=int(s[2])+len(s[6][:lo].replace('-',''));out.write(f's {s[1]} {newstart} {size} {s[4]} {s[5]} {text}\n')
   out.write('\n')
  with (Z/files[g]).open() as inp:
   for l in inp:
    if not l.strip():emit(block);block=[]
    else:block.append(l)
   emit(block)
 assert covered==1000,(g,covered)
 for name,model in models.items():
  wig=O/(g+'.'+name+'.wig');cmd=['/home/stardrako/miniforge3/envs/cactus/bin/phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(target)]
  with wig.open('w') as out,(O/(g+'.'+name+'.log')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=180)
  vals={};pos=None;step=1
  for l in wig.read_text().splitlines():
   if l.startswith('fixedStep'):
    fields=dict(x.split('=') for x in l.split()[1:]);pos=int(fields['start'])-1;step=int(fields.get('step',1))
   elif l and not l.startswith('#'):
    vals[pos]=float(l);pos+=step
  assert len(vals)==1000,(g,len(vals));v=np.array([vals[x] for x in sorted(vals)]);assert np.isfinite(v).all()
  rows.append(dict(genes=w['genes'],model=name,scored_bp=len(v),mean=float(v.mean()),max50bp_mean=float(np.convolve(v,np.ones(50)/50,'valid').max()),fraction_phyloP_ge_2=float(np.mean(v>=2))))
  print(rows[-1],flush=True)
with (O/'local_conservation_sensitivity.tsv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
trees=[re.search(r'TREE:\s*(.*)',p.read_text()).group(1) for p in models.values()];lengths=[np.array([float(x) for x in re.findall(r':([0-9.eE+-]+)',t)]) for t in trees];ratios=lengths[1]/lengths[0]
(O/'conservation_sensitivity_summary.json').write_text(json.dumps(dict(branch_scale_median=float(np.median(ratios)),branch_scale_min=float(ratios.min()),branch_scale_max=float(ratios.max()),results=rows,limitation='Second unfiltered random region scale sensitivity; not a fully independent neutral model or production-grade conservation validation.'),indent=2))
