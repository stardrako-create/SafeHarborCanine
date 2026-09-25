from pathlib import Path
import re,json,subprocess,time,numpy as np,csv
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.9';Z=Path('/mnt/d/Jin2024_work/zoonomia_hal');model=O/'Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod'
def leaves(p):return set(re.findall(r'[(,]([^():,]+):',re.search(r'TREE:\s*(.*)',p.read_text()).group(1)))
new=leaves(model);old=leaves(Z/'neutral_model_1region.mod');assert 'Canis_lupus_familiaris' in new
compat=dict(pilot_leaf_count=len(old),published_leaf_count=len(new),only_pilot=sorted(old-new),only_published=sorted(new-old));(O/'model_compatibility.json').write_text(json.dumps(compat,indent=2));print(compat,flush=True)
files={'ANO2':'maf_dog/cand_36.chr27.39320764.maf','LOC119876429':'new3_phylop/loc1_LOC119876429_LOC119872513.maf','NPNT':'new3_phylop/loc3_NPNT_TBCK.maf'}
results=[]
for gene,source in files.items():
 out=O/(gene+'.published_model.wig');cmd=['/home/stardrako/miniforge3/envs/cactus/bin/phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(Z/source)];t=time.time()
 with out.open('w') as f,(O/(gene+'.published_model.log')).open('w') as err:r=subprocess.run(cmd,stdout=f,stderr=err,timeout=600)
 assert r.returncode==0
 results.append(dict(gene=gene,command=cmd,seconds=time.time()-t,returncode=r.returncode));(O/'conservation_execution.json').write_text(json.dumps(results,indent=2));print(gene,'done',round(time.time()-t),flush=True)
