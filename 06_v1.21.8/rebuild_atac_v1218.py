from pathlib import Path
import numpy as np,pyBigWig,csv,json,time,math,sys
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.8';D=O/'ATAC';D.mkdir(exist_ok=True);(D/'cache').mkdir(exist_ok=True)
wrows=list(csv.DictReader((P/'06_v1.20.0/ATAC/qc_weights_unified76.tsv').open(),delimiter='\t'));samples=[r['sample'] for r in wrows];weights=np.array([float(r['weight']) for r in wrows]);handles=[];raws=[];scales=[];lambdas=[];events=[];t0=time.time()
for s in samples:
 root=Path('/mnt/d/Jin2024_work/ehsan5_pipeline/per_dog' if s.startswith('SRR307') else '/mnt/d/Jin2024_work/ATAC_pipeline/per_dog')/s
 qc=next(csv.DictReader((root/'qc'/f'{s}.qc.tsv').open(),delimiter='\t'));scales.append(float(qc['total_reads'])/1e6);lambdas.append(float(qc['lambda']));handles.append(pyBigWig.open(str(root/'bigwig'/f'{s}.cpm.bw')));rp=root/'bigwig'/f'{s}.raw.bw';raws.append(pyBigWig.open(str(rp)) if rp.exists() else None)
scales=np.array(scales);lambdas=np.array(lambdas);chroms=list(handles[0].chroms().items());assert all(h.chroms()==dict(chroms) for h in handles)
def bins(h,c,a,b):
 arr=np.zeros(math.ceil((b-a)/25)+1,dtype=np.float64);iv=h.intervals(c,a,b)
 if iv:
  x=np.array(iv,dtype=np.float64);s=np.maximum(x[:,0].astype(np.int64),a);e=np.minimum(x[:,1].astype(np.int64),b);v=x[:,2];assert np.isfinite(v).all() and (v>=0).all()
  assert np.all(s%25==0) and np.all((e%25==0)|(e==dict(chroms)[c])),(c,a,b,'unaligned intervals')
  np.add.at(arr,(s-a)//25,v);np.add.at(arr,(e-a+24)//25,-v)
 return np.cumsum(arr[:-1]).astype(np.float32)
def read(i,c,a,b):
 try:return bins(handles[i],c,a,b)
 except RuntimeError:
  try:return bins(handles[i],c,a,b)
  except RuntimeError as e:
   if raws[i] is None:raise RuntimeError(f'{samples[i]}:{c}:{a}-{b} no recovery source') from e
   v=bins(raws[i],c,a,b)/scales[i];events.append(dict(sample=samples[i],chrom=c,start=a,end=b,recovery='raw_bigwig_divided_by_total_reads_per_million'));(D/'recovery_events.json').write_text(json.dumps(events,indent=2));return v
if '--benchmark' in sys.argv:
 c='NC_051811.1';a=17235000;b=18235000;t=time.time();matrix=np.array([read(i,c,a,b) for i in range(len(samples))]);original=np.nan_to_num(handles[0].values(c,a,b,numpy=True),nan=0).reshape(-1,25).mean(1);assert np.allclose(matrix[0],original,rtol=1e-6,atol=1e-7)
 print(json.dumps(dict(chunk_seconds=time.time()-t,estimated_read_seconds=(sum(n for _,n in chroms)/1e6)*(time.time()-t),shape=matrix.shape,events=events)),flush=True);sys.exit()
try:
 for c,size in chroms:
  dest=D/'cache'/(c+'.npz')
  if dest.exists():continue
  n=math.ceil(size/25);means=np.zeros(n,dtype=np.float32);var=np.full(n,np.nan,dtype=np.float32);evidence=np.zeros(n,dtype=np.uint8)
  for a in range(0,size,2000000):
   b=min(size,a+2000000);mat=np.array([read(i,c,a,b) for i in range(len(samples))]);raw=mat*scales[:,None];gate=raw/(raw+lambdas[:,None])>=.2;wg=gate*weights[:,None];den=wg.sum(0);count=gate.sum(0);valid=den>1e-9;m=np.zeros(mat.shape[1],dtype=np.float32);m[valid]=((mat*wg).sum(0)[valid]/den[valid]).astype(np.float32)
   sort=np.ascontiguousarray(np.where(gate,mat,np.inf).T);sort.sort(axis=1);q=[]
   for p in [.25,.75]:
    ix=np.maximum(count-1,0)*p;l=np.floor(ix).astype(int);u=np.ceil(ix).astype(int);rr=np.arange(len(l));q.append(sort[rr,l]+(sort[rr,u]-sort[rr,l])*(ix-l))
   with np.errstate(invalid='ignore'):v=q[1]-q[0]
   v[count<2]=np.nan;z=a//25;means[z:z+len(m)]=m;var[z:z+len(m)]=v;evidence[z:z+len(m)]=count
  np.savez_compressed(dest,mean=means,variability=var,evidence=evidence);print(c,size,'elapsed',round(time.time()-t0),'recovered',len(events),flush=True)
 allmeans=[]
 for c,_ in chroms:
  with np.load(D/'cache'/(c+'.npz')) as d:allmeans.append(d['mean'][d['evidence']>0])
 background=float(np.median(np.concatenate(allmeans)));assert background>0;del allmeans
 output={k:pyBigWig.open(str(D/(n+'.partial')),'w') for k,n in [('mean','mother_track_accessibility_level.bw'),('variability','variability.bw'),('evidence','evidence_dogs.bw')]}
 for h in output.values():h.addHeader(chroms)
 for c,size in chroms:
  with np.load(D/'cache'/(c+'.npz')) as d:
   starts=np.arange(0,size,25);ends=np.minimum(starts+25,size)
   for k,h in output.items():
    v=d[k]/background if k=='mean' else d[k];mask=np.isfinite(v)
    if k=='mean':mask &= d['evidence']>0
    if mask.any():h.addEntries([c]*int(mask.sum()),starts[mask].tolist(),ends=ends[mask].tolist(),values=v[mask].astype(float).tolist())
 for h in output.values():h.close()
 for p in D.glob('*.bw.partial'):p.replace(p.with_suffix(''))
 (D/'build_manifest.json').write_text(json.dumps(dict(status='complete',samples=samples,chromosomes=len(chroms),genome_bp=sum(n for _,n in chroms),gain_background=background,elapsed_seconds=time.time()-t0,recoveries=events,mean_missing_policy='no gated dog means gap, not measured zero',variability_missing_policy='fewer than two gated dogs means gap',bin_size=25),indent=2))
except Exception as e:
 (D/'build_failure.json').write_text(json.dumps(dict(error=repr(e),elapsed_seconds=time.time()-t0,recoveries=events),indent=2));raise
