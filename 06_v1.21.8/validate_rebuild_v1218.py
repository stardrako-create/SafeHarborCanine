from pathlib import Path
import json,csv,warnings,math
import numpy as np,pyBigWig
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.8';D=O/'ATAC';manifest=json.loads((D/'build_manifest.json').read_text());gain=manifest['gain_background']
wrows=list(csv.DictReader((P/'06_v1.20.0/ATAC/qc_weights_unified76.tsv').open(),delimiter='\t'));weights=np.array([float(r['weight']) for r in wrows]);handles=[];rawhandles=[];scales=[];lambdas=[]
for row in wrows:
 s=row['sample'];root=Path('/mnt/d/Jin2024_work/ehsan5_pipeline/per_dog' if s.startswith('SRR307') else '/mnt/d/Jin2024_work/ATAC_pipeline/per_dog')/s;qc=next(csv.DictReader((root/'qc'/f'{s}.qc.tsv').open(),delimiter='\t'));scales.append(float(qc['total_reads'])/1e6);lambdas.append(float(qc['lambda']));handles.append(pyBigWig.open(str(root/'bigwig'/f'{s}.cpm.bw')));rp=root/'bigwig'/f'{s}.raw.bw';rawhandles.append(pyBigWig.open(str(rp)) if rp.exists() else None)
scales=np.array(scales);lambdas=np.array(lambdas);out={k:pyBigWig.open(str(D/f)) for k,f in [('mean','mother_track_accessibility_level.bw'),('var','variability.bw'),('evidence','evidence_dogs.bw')]};chroms=handles[0].chroms();assert all(h.chroms()==chroms for h in out.values())
def dense(h,c,a,b):
 v=np.nan_to_num(h.values(c,a,b,numpy=True),nan=0).astype(np.float32);n=math.ceil(len(v)/25);pad=n*25-len(v);v=np.pad(v,(0,pad));s=v.reshape(-1,25).sum(1,dtype=np.float64);length=np.minimum(25,np.arange(n)*-25+(b-a));return (s/length).astype(np.float32)
recovery=[]
for event in manifest['recoveries']:
 i=[r['sample'] for r in wrows].index(event['sample']);c=event['chrom'];matches=[];failures=[]
 for a in np.linspace(event['start'],event['end']-1000,9).astype(int):
  a=int(a//25*25)
  try:
   cp=dense(handles[i],c,a,a+1000);raw=dense(rawhandles[i],c,a,a+1000)/scales[i];assert np.allclose(cp,raw,atol=2e-6,rtol=2e-5),(event,a,float(np.max(abs(cp-raw))))
   matches.append(dict(start=a,max_abs_error=float(np.max(abs(cp-raw)))))
  except RuntimeError:failures.append(a)
 assert len(matches)>=3,(event,'insufficient readable comparison')
 recovery.append(dict(event=event,readable_checks=matches,unreadable_checks=failures))
regions=[]
for row in csv.DictReader((P/'06_v1.21.5/cohort_and_support.tsv').open(),delimiter='\t'):
 if row['genes'].endswith('_region'):continue
 a=int(row['start'])//25*25;regions.append((row['chrom'],a,a+1025,'selected_'+row['genes']))
for e in manifest['recoveries']:
 a=e['start']+1000000;regions.append((e['chrom'],a,a+1000,'recovered_region'))
rng=np.random.default_rng(1218);main=[(c,n) for c,n in chroms.items() if c.startswith('NC_') and n>2000000]
for _ in range(12):
 c,size=main[int(rng.integers(len(main)))];a=int(rng.integers(size-2000))//25*25;regions.append((c,a,a+1000,'random'))
c,size=next((c,n) for c,n in chroms.items() if n%25 and n<2000);regions.append((c,0,size,'partial_final_bin'))
checks=[]
for c,a,b,label in regions:
 mat=[]
 for i in range(len(handles)):
  try:v=dense(handles[i],c,a,b)
  except RuntimeError:v=dense(rawhandles[i],c,a,b)/scales[i]
  mat.append(v)
 mat=np.array(mat);raw=mat*scales[:,None];gate=raw/(raw+lambdas[:,None])>=.2;count=gate.sum(0);wg=gate*weights[:,None];den=wg.sum(0);mean=np.full(len(den),np.nan);valid=den>1e-9;mean[valid]=(mat*wg).sum(0)[valid]/den[valid]/gain
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',RuntimeWarning);q=np.nanpercentile(np.where(gate,mat,np.nan),[25,75],axis=0);var=q[1]-q[0]
 var[count<2]=np.nan
 for k,expected in [('mean',mean),('var',var),('evidence',count)]:
  observed=out[k].values(c,a,b,numpy=True)[::25];assert np.array_equal(np.isnan(observed),np.isnan(expected)),(c,a,k,'missingness');assert np.allclose(observed,expected,rtol=2e-5,atol=2e-6,equal_nan=True),(c,a,k,float(np.nanmax(abs(observed-expected))))
 checks.append(dict(chrom=c,start=a,end=b,label=label,bins=len(count),dense_mean_IQR_evidence_match=True))
result=dict(status='passed',recovery_scale_checks=recovery,independent_dense_checks=checks,total_independent_bins=sum(x['bins'] for x in checks),limitations='Targeted independent numerical checks, not read-level validation of source BigWigs or absolute ATAC QC. Whole-track cache invariants checked separately.')
(O/'rebuild_targeted_validation.json').write_text(json.dumps(result,indent=2));print('PASS',len(checks),'regions',result['total_independent_bins'],'bins; four recovery scale checks')
