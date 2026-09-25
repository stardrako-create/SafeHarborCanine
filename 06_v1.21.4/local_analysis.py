from pathlib import Path
import csv,json,hashlib,shutil,sys,gzip,collections
import numpy as np
import pyBigWig
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
Z=Path('/mnt/d/Jin2024_work/zoonomia_hal')
O=P/'06_v1.21.4'; O.mkdir(exist_ok=True)
def table(path,rows):
 if not rows:return
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
rows=read(P/'06_v1.21.3/candidates_scored_v1213.tsv')
candidates=[r for r in rows if r['hard_veto']=='False']
for r in candidates:r['start']=int(r['start']);r['end']=int(r['end'])
shutil.copy2(__file__,O/'local_analysis.py')
sys.path.insert(0,str(P/'scripts'))
from bw_utils import evidence_summary
manifest=json.loads((P/'06_v1.21.1/manifest.json').read_text())
args=manifest['command']; paths={args[i][2:]:Path(args[i+1]) for i in range(2,len(args)-1,2)}
tracks={k:pyBigWig.open(str(paths[k])) for k in ['atac-mean-bw','atac-variability-bw','atac-peak-frequency-bw','rrbs-mean-bw','rrbs-coverage-bw']}
def intervals(path):
 d=collections.defaultdict(list)
 for line in path.open():
  if line.startswith('#') or not line.strip():continue
  a=line.split()
  if a[1:3]==['start','end']:continue
  d[a[0]].append((int(a[1]),int(a[2])))
 return d
annotations={k:intervals(paths[k]) for k in ['atac-peaks-bed','external-regulatory-bed','tad-boundaries-bed']}
annotations['cpg']=intervals(P/'06_v1.21.1/cpg_islands_independent_calls.bed')
repeats=collections.defaultdict(list)
for rp in [Z/'new3_repeatmasker/new3_candidates.fa.out',Path('/mnt/d/Jin2024_work/repeatmasker_43/candidates_43.fa.out')]:
 for line in rp.open():
  a=line.split()
  if not a or not a[0].isdigit():continue
  for c in candidates:
   query=f"{c['chrom']}:{c['start']}-{c['end']}"
   if a[4]==query or a[4].endswith('_'+c['left_gene']+'_'+c['right_gene']):
    repeats[c['left_gene']].append((c['start']+int(a[5])-1,c['start']+int(a[6])))
def union_length(ivals):
 end=-1;total=0
 for s,e in sorted(ivals):total+=max(0,e-max(s,end));end=max(end,e)
 return total
def overlap_bp(ivals,s,e):return union_length([(max(a,s),min(b,e)) for a,b in ivals if a<e and b>s])
# Chain evidence is alignment coverage, not proof of biological rearrangement.
chain_path=Z/'liftover/GCF_014441545.1ToCanFam3.over.chain.gz'
chainhits=collections.defaultdict(list);h=None;blocks=[]
def finish_chain():
 if h is None:return
 for c in candidates:
  if c['chrom']!=h[2]:continue
  selected=[]
  for ts,te,qs,qe in blocks:
   a=max(ts,c['start']);b=min(te,c['end'])
   if a<b:selected.append((a,b,qs+(a-ts),qs+(b-ts)))
  if selected:chainhits[c['chrom'],c['start']].append({'id':h[12],'score':int(h[1]),'target':h[7],'strand':h[9],'qsize':int(h[8]),'covered':sum(b-a for a,b,_,_ in selected),'blocks':selected})
with gzip.open(chain_path,'rt') as f:
 for line in f:
  a=line.split()
  if not a:continue
  if a[0]=='chain':
   finish_chain();h=a;tp=int(a[5]);qp=int(a[10]);blocks=[]
  else:
   n=int(a[0]);blocks.append((tp,tp+n,qp,qp+n));tp+=n;qp+=n
   if len(a)==3:tp+=int(a[1]);qp+=int(a[2])
 finish_chain()
chain_report=[]
for c in candidates:
 hits=sorted(chainhits[c['chrom'],c['start']],key=lambda x:x['covered'],reverse=True)
 c['_chain']=hits[0] if hits else None
 chain_report.append({'genes':c['left_gene']+'/'+c['right_gene'],'region_bp':c['end']-c['start'],'chains_with_aligned_bases':len(hits),'best_chain_covered_bp':hits[0]['covered'] if hits else 0,'best_chain_coverage_fraction':hits[0]['covered']/(c['end']-c['start']) if hits else 0,'best_chain_target':hits[0]['target'] if hits else '', 'best_chain_strand':hits[0]['strand'] if hits else '', 'targets_any_chain':';'.join(sorted({x['target'] for x in hits}))})
table(O/'chain_alignment_audit.tsv',chain_report)
(O/'chain_blocks.json').write_text(json.dumps({c['left_gene']:chainhits[c['chrom'],c['start']] for c in candidates},indent=2))
print('Chain audit complete',flush=True)
def wig_values(p):
 values={};chrom=None;pos=0;step=1
 for line in p.open():
  if line.startswith('fixedStep'):
   h=dict(x.split('=') for x in line.split()[1:]);chrom=h['chrom'];pos=int(h['start'])-1;step=int(h.get('step',1))
  elif line.startswith(('track','#')) or not line.strip():continue
  else:values[pos]=float(line.strip());pos+=step
 return values
mapped={};cons=[]
for c in candidates:
 gene=c['left_gene'];wig=None
 if gene=='LOC119876429':wig=Z/'new3_phylop/loc1_LOC119876429_LOC119872513.wig'
 if gene=='NPNT':wig=Z/'new3_phylop/loc3_NPNT_TBCK.wig'
 if gene=='ANO2':
  files=list((Z/'phylop_out').glob('*.chr27.39320764.wig'));wig=files[0] if len(files)==1 else None
 if wig is None:continue
 vals=wig_values(wig);chain=c['_chain'];arr=np.full(c['end']-c['start'],np.nan)
 for a,b,qs,qe in chain['blocks']:
  for t,q in zip(range(a,b),range(qs,qe)):
   q=q if chain['strand']=='+' else chain['qsize']-1-q
   if q in vals:arr[t-c['start']]=vals[q]
 mapped[gene]=arr
 v=np.array(list(vals.values()));valid=np.isfinite(v)
 roll=np.convolve(np.where(valid,v,0),np.ones(50)/50,'valid');support=np.convolve(valid.astype(int),np.ones(50,dtype=int),'valid')
 maximum=float(np.max(roll[support==50])) if np.any(support==50) else None
 cons.append({'genes':gene+'/'+c['right_gene'],'wig':str(wig),'wig_values':len(v),'recomputed_max50':maximum,'recorded_max50':c['max_50bp_rolling_phyloP'],'max50_matches':abs(maximum-float(c['max_50bp_rolling_phyloP']))<1e-6,'source_region_scored_fraction':float(np.isfinite(arr).mean()),'model':'pilot_1region_100kb_not_revalidated'})
table(O/'conservation_audit.tsv',cons)
# Fixed local resolutions; all windows retained. Background measures matching length.
rng=np.random.default_rng(42);bw=tracks['atac-mean-bw'];chroms=bw.chroms();names=[k for k,v in chroms.items() if v>=2000 and k.startswith('NC_')];sizes=np.array([chroms[k] for k in names]);bg={}
for width in (500,1000,2000):
 samples=[];weights=sizes-width+1;weights=weights/weights.sum()
 for _ in range(3000):
  chrom=str(rng.choice(names,p=weights));start=int(rng.integers(0,chroms[chrom]-width+1));v=bw.stats(chrom,start,start+width,type='mean',exact=True)[0]
  if v is not None and np.isfinite(v):samples.append(v)
 bg[width]=np.sort(samples)
table(O/'local_backgrounds.tsv',[{'window_bp':k,'observed_n':len(v),'p50':float(np.percentile(v,50)),'p55':float(np.percentile(v,55)),'p90':float(np.percentile(v,90))} for k,v in bg.items()])
local=[]
for c in candidates:
 if c['evaluation_status']!='passes_recorded_checks':continue
 ch,s,e=c['chrom'],c['start'],c['end'];gene=c['left_gene'];arrs={k:np.array(b.values(ch,s,e)) for k,b in tracks.items()}
 for width in (500,1000,2000):
  starts=sorted(set(range(s,e-width+1,250))|{e-width})
  for a in starts:
   b=a+width;sl=slice(a-s,b-s);v=arrs['atac-mean-bw'][sl];finite=np.isfinite(v);mean=float(v[finite].mean()) if finite.any() else None
   meth=arrs['rrbs-mean-bw'][sl];cov=arrs['rrbs-coverage-bw'][sl];mask=np.isfinite(meth)&np.isfinite(cov)&(cov>=1)
   phy=mapped.get(gene,np.full(e-s,np.nan))[sl]
   freq=arrs['atac-peak-frequency-bw'][sl];fv=freq[np.isfinite(freq)]
   r={'genes':gene+'/'+c['right_gene'],'chrom':ch,'start':a,'end':b,'window_bp':width,'atac_mean':mean,'atac_observed_fraction':float(finite.mean()),'atac_length_matched_percentile':float(np.searchsorted(bg[width],mean,side='right')/len(bg[width])*100) if mean is not None else None,'atac_peak_frequency_max':float(fv.max()) if len(fv) else 0.0,'rrbs_observed_fraction':float(mask.mean()),'rrbs_mean_observed':float(meth[mask].mean()) if mask.any() else None,'rrbs_dogs_max':float(np.nanmax(cov)) if np.isfinite(cov).any() else None,'repeat_overlap_bp':overlap_bp(repeats[gene],a,b) if repeats[gene] else None,'phylop_scored_fraction':float(np.isfinite(phy).mean()),'phylop_mean_observed':float(np.nanmean(phy)) if np.isfinite(phy).any() else None,'phylop_max_observed':float(np.nanmax(phy)) if np.isfinite(phy).any() else None}
   for name,ann in annotations.items():r[name+'_overlap_bp']=overlap_bp(ann.get(ch,[]),a,b)
   local.append(r)
table(O/'all_local_windows.tsv',local)
tops=[]
for gene in sorted({r['genes'] for r in local}):
 for width in (500,1000,2000):
  group=[r for r in local if r['genes']==gene and r['window_bp']==width]
  chosen=[]
  for r in sorted(group,key=lambda x:x['atac_mean'] if x['atac_mean'] is not None else -1,reverse=True):
   if all(r['end']<=q['start'] or r['start']>=q['end'] for q in chosen):chosen.append(r)
   if len(chosen)==3:break
  tops+=chosen
table(O/'top_local_windows_exploratory.tsv',tops)
used=[P/'06_v1.21.3/candidates_scored_v1213.tsv',chain_path,Z/'neutral_model_1region.mod',P/'06_v1.21.1/cpg_islands_independent_calls.bed']
(O/'analysis_manifest.json').write_text(json.dumps({'version':'1.21.4','seed':42,'grid_step_bp':250,'widths_bp':[500,1000,2000],'background_attempts_per_width':3000,'background_scope':'NC_ reference chromosomes; length matched, not other covariates','selection':'top local windows by ATAC are exploratory, not safe insertion nominations','track_paths':{k:str(paths[k]) for k in tracks},'hashes':{str(p):sha(p) for p in used}},indent=2))
print(json.dumps({'local_windows':len(local),'top_windows':len(tops),'conservation':cons,'chain':chain_report},indent=2),flush=True)
