from pathlib import Path
import numpy as np,json,csv
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.9';Z=Path('/mnt/d/Jin2024_work/zoonomia_hal');chains=json.loads((P/'06_v1.21.4/chain_blocks.json').read_text());oldfiles={'ANO2':Z/'phylop_out/cand_36.chr27.39320764.wig','LOC119876429':Z/'new3_phylop/loc1_LOC119876429_LOC119872513.wig','NPNT':Z/'new3_phylop/loc3_NPNT_TBCK.wig'}
def readwig(path):
 vals={}
 for line in path.open():
  if line.startswith('fixedStep'):
   h=dict(x.split('=') for x in line.split()[1:]);pos=int(h['start'])-1;step=int(h.get('step',1))
  elif line.strip() and not line.startswith(('#','track')):vals[pos]=float(line);pos+=step
 return vals
regions={}
for r in csv.DictReader((P/'06_v1.21.8/candidates_scored_v1218.tsv').open(),delimiter='\t'):
 if r['left_gene'] in oldfiles:regions[r['left_gene']]=r
arrays={};summaries=[]
def summarize(v):
 finite=np.isfinite(v);ok=np.convolve(finite.astype(int),np.ones(50,dtype=int),'valid')==50;means=np.convolve(np.nan_to_num(v,nan=0),np.ones(50)/50,'valid');return dict(bp=len(v),observed_fraction=float(finite.mean()),mean=float(v[finite].mean()),max50bp_mean=float(means[ok].max()) if ok.any() else None,fraction_observed_ge2=float(np.mean(v[finite]>=2)))
for gene,r in regions.items():
 a=int(r['start']);b=int(r['end']);c=max(chains[gene],key=lambda x:x['covered']);assert c['strand']=='+';new=readwig(O/(gene+'.published_model.wig'));old=readwig(oldfiles[gene]);assert set(new)==set(old),(gene,'coordinate mismatch')
 arrs={}
 for label,d in [('pilot',old),('published',new)]:
  v=np.full(b-a,np.nan)
  for ts,te,qs,qe in c['blocks']:
   for t,q in zip(range(ts,te),range(qs,qe)):
    if a<=t<b and q in d:v[t-a]=d[q]
  arrs[label]=v
  summaries.append(dict(gene=gene,scope='region',model=label,**summarize(v)))
 arrays[gene]=(a,arrs)
 np.savez_compressed(O/(gene+'.conservation_ROS.npz'),start=a,**arrs)
local=[]
for r in csv.DictReader((P/'06_v1.21.8/local_windows_evidence_status.tsv').open(),delimiter='\t'):
 gene=r['genes'].split('/')[0];a,arrs=arrays[gene];s=int(r['start'])-a;e=int(r['end'])-a
 for label,v in arrs.items():
  stats=summarize(v[s:e]);r[label+'_phyloP_observed_fraction']=stats['observed_fraction'];r[label+'_phyloP_max50bp_mean']=stats['max50bp_mean'];r[label+'_phyloP_mean']=stats['mean']
 r['conservation_status']='published_ancestral_repeat_model_comparison_available_not_safe_harbor_validation';local.append(r)
 if r['selected_1kb_previously']=='True':
  for label,v in arrs.items():summaries.append(dict(gene=gene,scope='previous_selected_1kb',model=label,**summarize(v[s:e])))
def write(path,rows):
 with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write(O/'conservation_model_comparison.tsv',summaries);write(O/'local_windows_conservation_rechecked.tsv',local)
(O/'conservation_review.json').write_text(json.dumps(dict(model_source='https://cgl.gi.ucsc.edu/data/cactus/241-mammals-human-phylop-models/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod',model_training='General model from ancestral repeat positions, per source README',matching_leaf_sets=241,local_windows=len(local),summary=summaries,limitations=['Same underlying alignment; independent neutral calibration, not independent biological data','No genome-wide FDR calculation reproduced','Old max50bp threshold 6.5 not validated under this model; no automatic safety pass','Only three full regions scored with published model in this run','No new functional or T-cell validation']),indent=2));print(json.dumps(summaries,indent=2))
