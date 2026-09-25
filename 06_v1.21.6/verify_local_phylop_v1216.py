from pathlib import Path
import csv,json,numpy as np
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.6';Z=Path('/mnt/d/Jin2024_work/zoonomia_hal')
paths={'ANO2':'phylop_out/cand_36.chr27.39320764.wig','LOC119876429':'new3_phylop/loc1_LOC119876429_LOC119872513.wig','NPNT':'new3_phylop/loc3_NPNT_TBCK.wig'}
def load(path):
 vals={};pos=None
 for l in path.open():
  if l.startswith('fixedStep'):
   d=dict(x.split('=') for x in l.split()[1:]);pos=int(d['start'])-1;step=int(d.get('step',1))
  elif l.strip() and not l.startswith(('#','track','variableStep')):
   vals[pos]=float(l);pos+=step
 return vals
out=[]
for r in csv.DictReader((P/'06_v1.21.5/local_canfam3_blocks.tsv').open(),delimiter='\t'):
 g=r['genes'].split('/')[0];old=load(Z/paths[g]);new=load(O/(g+'.pilot.wig'));a=int(r['start']);b=int(r['end'])
 assert set(new)==set(range(a,b)),(g,min(new),max(new),a,b)
 diffs=[abs(new[x]-old[x]) for x in range(a,b)];out.append(dict(genes=r['genes'],bases_compared=len(diffs),max_abs_difference=max(diffs),identical_at_reported_precision=all(d==0 for d in diffs)))
(O/'local_phylop_reproduction.json').write_text(json.dumps(out,indent=2));print(out)
