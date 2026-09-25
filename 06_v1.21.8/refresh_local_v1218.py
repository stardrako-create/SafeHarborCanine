from pathlib import Path
import csv,json,collections
import numpy as np,pyBigWig
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.8';D=O/'ATAC'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
 with p.open('w') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
rows=read(P/'06_v1.21.7/local_windows_evidence_status.tsv');bw=pyBigWig.open(str(D/'mother_track_accessibility_level.bw'));ev=pyBigWig.open(str(D/'evidence_dogs.bw'));chroms=bw.chroms();names=[k for k,v in chroms.items() if v>=2000 and k.startswith('NC_')];sizes=np.array([chroms[k] for k in names]);rng=np.random.default_rng(42);bg={};bgrows=[]
for width in [500,1000,2000]:
 values=[];prob=sizes-width+1;prob=prob/prob.sum()
 for i in range(3000):
  c=str(rng.choice(names,p=prob));a=int(rng.integers(0,chroms[c]-width+1));v=bw.stats(c,a,a+width,type='mean',exact=True)[0];bgrows.append(dict(width=width,chrom=c,start=a,end=a+width,value=v))
  if v is not None and np.isfinite(v):values.append(v)
 bg[width]=np.sort(values)
write(O/'local_background_windows.tsv',bgrows)
for r in rows:
 c=r['chrom'];a=int(r['start']);b=int(r['end']);v=bw.values(c,a,b,numpy=True);finite=np.isfinite(v);n=ev.values(c,a,b,numpy=True);mean=float(np.mean(v[finite],dtype=np.float64)) if finite.any() else None
 r['previous_atac_mean']=r['atac_mean'];r['previous_atac_length_matched_percentile']=r['atac_length_matched_percentile'];r['atac_mean']=mean;r['atac_observed_fraction']=float(finite.mean());r['atac_length_matched_percentile']=float(np.searchsorted(bg[b-a],mean,side='right')/len(bg[b-a])*100) if mean is not None else None;r['gate_evidence_dogs_min']=int(n.min());r['gate_evidence_dogs_median']=float(np.median(n));r['gate_evidence_dogs_max']=int(n.max());r['exploratory_Pareto_front']=False
# Retain explicit trade-offs rather than choosing new score weights or a convenient ATAC cutoff.
front=[];eligible_counts={}
for gene in sorted({r['genes'] for r in rows}):
 group=[r for r in rows if r['genes']==gene and int(r['window_bp'])==1000 and int(r['epic_promoter_enhancer_bp_union'])==0 and int(r['external-regulatory-bed_overlap_bp'])==0 and r['atac_mean'] is not None and all(r[k] not in ['',None] for k in ['repeat_overlap_bp','rrbs_observed_fraction','rrbs_mean_observed'])];eligible_counts[gene]=len(group)
 x=np.array([[float(r['atac_mean']),-float(r['repeat_overlap_bp']),float(r['rrbs_observed_fraction']),-float(r['rrbs_mean_observed'])] for r in group])
 for i,r in enumerate(group):
  dominated=np.any(np.all(x>=x[i],axis=1)&np.any(x>x[i],axis=1))
  if not dominated:r['exploratory_Pareto_front']=True;front.append(r)
write(O/'local_windows_evidence_status.tsv',rows)
if front:write(O/'local_tradeoff_frontier_exploratory.tsv',front)
selected=[{k:r[k] for k in ['genes','start','end','atac_mean','atac_length_matched_percentile','previous_atac_length_matched_percentile','atac_observed_fraction','gate_evidence_dogs_min','gate_evidence_dogs_median','gate_evidence_dogs_max','epic_promoter_enhancer_bp_union','decision_status']} for r in rows if r['selected_1kb_previously']=='True']
(O/'local_refresh_summary.json').write_text(json.dumps(dict(local_windows=len(rows),background_observed_n={str(k):len(v) for k,v in bg.items()},previously_selected=selected,eligible_1kb_without_recorded_regulatory_overlap=eligible_counts,frontier_counts=dict(collections.Counter(r['genes'] for r in front)),frontier_definition='Within each previously selected region: maximize ATAC measured mean and RRBS observed fraction; minimize repeat bp and observed methylation. No new safety threshold or validation implied.',limitations=['Only 3 regions examined','Missing covariate windows omitted from trade-off comparison, not biologically excluded','Adjacent/overlapping windows are not independent','No correction for exploratory search of local maxima','Most alternatives lack local variant, remapping and T-context checks','Absence of mapped regulatory annotation is not proof of neutrality']),indent=2));print(json.dumps(selected,indent=2));print('Frontier',len(front),eligible_counts)
