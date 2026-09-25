from pathlib import Path
import csv,json,hashlib,sys,shutil
import numpy as np
import pyBigWig
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.5';O.mkdir(exist_ok=True)
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
weights=read(P/'06_v1.20.0/ATAC/qc_weights_unified76.tsv');assert len(weights)==76
local=read(P/'06_v1.21.4/all_local_windows.tsv')
windows=[]
for gene in sorted({x['genes'] for x in local}):
 r=max([x for x in local if x['genes']==gene and x['window_bp']=='1000'],key=lambda x:float(x['atac_mean']))
 windows.append({'genes':gene,'chrom':r['chrom'],'start':int(r['start']),'end':int(r['end'])})
regions=read(P/'06_v1.21.3/candidates_scored_v1213.tsv');regions=[r for r in regions if r['evaluation_status']=='passes_recorded_checks']
for r in regions:windows.append({'genes':r['left_gene']+'/'+r['right_gene']+'_region','chrom':r['chrom'],'start':int(r['start']),'end':int(r['end'])})
data={r['genes']:[] for r in windows};per=[];qc=[];errors=[];handles=[]
for idx,w in enumerate(weights):
 s=w['sample'];cohort='original71' if s.startswith('SRR273') else 'additional5'
 bases=[Path('/mnt/d/Jin2024_work/ATAC_pipeline/per_dog'),Path('/mnt/d/Jin2024_work/ehsan5_pipeline/per_dog')]
 dirs=[b/s for b in bases if (b/s/'bigwig'/f'{s}.cpm.bw').exists()];assert dirs,s
 d=dirs[0];q=read(d/'qc'/f'{s}.qc.tsv')[0];scale=float(q['total_reads'])/1e6;lam=float(q['lambda']);assert lam>0 and scale>0
 qc.append({'sample':s,'cohort':cohort,'weight':float(w['weight']),'frip':float(q['frip']),'tss_enrichment':float(q['tss_enrichment']),'reads':int(q['total_reads']),'lambda':lam,'path':str(d)})
 bw=pyBigWig.open(str(d/'bigwig'/f'{s}.cpm.bw'))
 peaks=[]
 for line in (d/'peaks'/f'{s}_peaks.narrowPeak').open():
  a=line.split();peaks.append((a[0],int(a[1]),int(a[2])))
 for r in windows:
  a=r['start']//25*25;b=((r['end']+24)//25)*25
  try:v=np.array(bw.values(r['chrom'],a,b,numpy=True),dtype=float)
  except RuntimeError as e:errors.append({'sample':s,'window':r['genes'],'error':str(e)});raise
  # CPM bigWig gaps mean zero coverage; unreadable data are errors, not zeros.
  v=np.nan_to_num(v,nan=0).reshape(-1,25).mean(axis=1);raw=v*scale;gate=raw/(raw+lam)>=.2
  lengths=np.array([max(0,min(t+25,r['end'])-max(t,r['start'])) for t in range(a,b,25)])
  cpm=float(np.average(v,weights=lengths));peak=any(ch==r['chrom'] and st<r['end'] and en>r['start'] for ch,st,en in peaks)
  per.append({'sample':s,'cohort':cohort,**r,'cpm_mean':cpm,'raw_to_genome_mean_ratio':cpm*scale/lam,'gate_fraction':float(np.average(gate,weights=lengths)),'has_peak_overlap':peak})
  data[r['genes']].append((v,gate,lengths))
 bw.close()
 print(f'{idx+1}/76 {s}',flush=True) if (idx+1)%10==0 else None
assert sum(q['cohort']=='original71' for q in qc)==71
write(O/'per_dog_local_atac.tsv',per);write(O/'qc_by_dog.tsv',qc)
w=np.array([q['weight'] for q in qc]);i71=np.array([q['cohort']=='original71' for q in qc]);summ=[]
def gated(mat,gates,weight,lengths):
 effective=gates*weight[:,None];den=effective.sum(axis=0);means=np.divide((mat*effective).sum(axis=0),den,out=np.full(mat.shape[1],np.nan),where=den>0)
 valid=np.isfinite(means);return float(np.average(means[valid],weights=lengths[valid])) if valid.any() else np.nan
for r in windows:
 name=r['genes'];v=np.stack([a[0] for a in data[name]]);g=np.stack([a[1] for a in data[name]]);lengths=data[name][0][2]
 m76=gated(v,g,w,lengths);m71=gated(v[i71],g[i71],w[i71],lengths)
 relevant=[x for x in per if x['genes']==name];ratios=np.array([x['raw_to_genome_mean_ratio'] for x in relevant])
 summ.append({**r,'gated_cpm_76':m76,'gated_cpm_71_fixed_qc_weights':m71,'additional5_change_pct':100*(m76/m71-1),'gate_dogs_min':int(g.sum(axis=0).min()),'gate_dogs_median':float(np.median(g.sum(axis=0))),'gate_dogs_max':int(g.sum(axis=0).max()),'dogs_with_peak_71':sum(x['has_peak_overlap'] and x['cohort']=='original71' for x in relevant),'dogs_with_peak_5':sum(x['has_peak_overlap'] and x['cohort']=='additional5' for x in relevant),'dogs_above_own_genome_mean':int((ratios>1).sum()),'raw_background_ratio_median':float(np.median(ratios)),'raw_background_ratio_q25':float(np.percentile(ratios,25)),'raw_background_ratio_q75':float(np.percentile(ratios,75))})
write(O/'cohort_and_support.tsv',summ)
# Paired bootstrap over original dogs, fixed previously selected windows.
rng=np.random.default_rng(42);ids=np.flatnonzero(i71);boot={r['genes']:[] for r in windows[:3]}
for j in range(2000):
 sampled=rng.choice(ids,size=71,replace=True)
 for r in windows[:3]:
  name=r['genes'];d=data[name];v=np.stack([d[i][0] for i in sampled]);g=np.stack([d[i][1] for i in sampled]);boot[name].append(gated(v,g,w[sampled],d[0][2]))
bootstrap=[]
for name,vals in boot.items():bootstrap.append({'genes':name,'bootstrap_n':2000,'gated_CPM_p025':float(np.percentile(vals,2.5)),'gated_CPM_p50':float(np.percentile(vals,50)),'gated_CPM_p975':float(np.percentile(vals,97.5))})
write(O/'fixed_window_bootstrap.tsv',bootstrap)
pair=[];names=list(boot)
for a in range(3):
 for b in range(a+1,3):
  diff=np.array(boot[names[a]])-np.array(boot[names[b]])
  pair.append({'a':names[a],'b':names[b],'difference_p025':float(np.percentile(diff,2.5)),'difference_p975':float(np.percentile(diff,97.5)),'fraction_a_above_b':float((diff>0).mean())})
write(O/'paired_bootstrap_differences.tsv',pair)
technical=[]
for r in windows[:3]:
 name=r['genes'];d=data[name];v=np.stack([a[0] for a in d]);lengths=d[0][2]
 scale=np.array([q['reads']/1e6 for q in qc]);lam=np.array([q['lambda'] for q in qc]);raw=v*scale[:,None];conf=raw/(raw+lam[:,None])
 np.savez_compressed(O/(name.split('/')[0]+'_per_dog_bins.npz'),cpm=v,confidence=conf,lengths=lengths,samples=np.array([q['sample'] for q in qc]))
 frip=np.array([q['frip'] for q in qc]);tss=np.array([q['tss_enrichment'] for q in qc])
 high=i71&(frip>=np.median(frip[i71]))&(tss>=np.median(tss[i71]))
 for label,subset,wt in [('71_QC_weighted',i71,w),('71_equal_weight',i71,np.ones(76)),('71_above_both_QC_medians',high,w)]:
  for floor in [.1,.2,.3]:
   value=gated(v[subset],conf[subset]>=floor,wt[subset],lengths)
   technical.append({'genes':name,'configuration':label,'included_dogs':int(subset.sum()),'confidence_floor':floor,'gated_CPM':value})
write(O/'technical_sensitivity.tsv',technical)
(O/'execution.json').write_text(json.dumps({'samples':76,'original':71,'additional':5,'windows':windows,'read_errors':errors,'bootstrap_seed':42,'note':'Fixed-window descriptive bootstrap; does not correct prior window selection; CPM before cohort-specific gain. QC weights held fixed when comparing 71/76.'},indent=2))
shutil.copy2(__file__,O/'clean_v1215.py')
print(json.dumps({'support':summ,'bootstrap':bootstrap,'pair':pair},indent=2),flush=True)
