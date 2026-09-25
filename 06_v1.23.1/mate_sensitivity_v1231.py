from pathlib import Path
import csv,json,gzip,subprocess,time,hashlib,shutil,fcntl,os,traceback
import numpy as np
from scipy.stats import spearmanr
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.23.1';O.mkdir(exist_ok=True)
lock=open('/tmp/safeharbor_v1231.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
rows=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'));idx=P/'06_v1.22.7/STAR_index';STAR='/home/stardrako/miniforge3/envs/cart_rna/bin/STAR';N=100000
state={'status':'starting','completed':[],'started':time.time()}
def write(p,x):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(x,indent=2));t.replace(p)
def save():state['updated']=time.time();write(O/'status.json',state)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
expected=[l.split('\t')[0] for l in (idx/'geneInfo.tab').read_text().splitlines()[1:]];assert len(expected)==42309
reference=json.loads((P/'06_v1.22.7/index_validated.json').read_text());assert sha(idx/'geneInfo.tab')==reference['index_sha256']['geneInfo.tab']
metrics=[];comparisons=[]
try:
 for r in rows:
  run=r['SRR'];folder=O/run;folder.mkdir(exist_ok=True);inputs=[folder/f'first100k_R{i}.fastq' for i in [1,2]]
  if shutil.disk_usage(O).free<20_000_000_000 or shutil.disk_usage('/mnt/c').free<15_000_000_000:raise RuntimeError('Insufficient disk reserve')
  state.update(status='extracting_subset',sample=run);save()
  # All three modes use exactly the same first 100,000 treated pairs, no additional filtering.
  with gzip.open(P/f'06_v1.22.6/samples/{run}/trimmed_1.fastq.gz','rt') as a,gzip.open(P/f'06_v1.22.6/samples/{run}/trimmed_2.fastq.gz','rt') as b,inputs[0].open('w') as x,inputs[1].open('w') as y:
   for k in range(N):
    ra=[a.readline() for _ in range(4)];rb=[b.readline() for _ in range(4)]
    assert ra[0].startswith('@') and rb[0].startswith('@') and ra[2].startswith('+') and rb[2].startswith('+')
    assert len(ra[1].strip())==len(ra[3].strip()) and len(rb[1].strip())==len(rb[3].strip())
    assert ra[0].split()[0].removesuffix('/1')==rb[0].split()[0].removesuffix('/2')
    x.writelines(ra);y.writelines(rb)
  hashes=[sha(p) for p in inputs];vectors={}
  for mode,reads in [('paired',inputs),('R1',[inputs[0]]),('R2',[inputs[1]])]:
   out=folder/mode;out.mkdir(exist_ok=True);state.update(status='aligning',mode=mode);save()
   cmd=[STAR,'--runThreadN','6','--genomeDir',str(idx),'--readFilesIn',*[str(p) for p in reads],'--outSAMtype','None','--quantMode','GeneCounts','--outFileNamePrefix',str(out/'STAR_'),'--outTmpDir',f'/tmp/sh_v1231_{run}_{mode}_{os.getpid()}']
   signature={'input_sha256':hashes,'mode':mode,'N':N,'reference_signature':reference['signature'],'command_without_temp':cmd[:-1]};receipt=out/'validated.json'
   if receipt.exists():assert json.loads(receipt.read_text())['signature']==signature
   else:
    start=time.time()
    with (out/'console.log').open('w') as f:ret=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    if ret.returncode:raise RuntimeError(f'{run}/{mode} STAR failed: {ret.returncode}')
   assert 'ALL DONE!' in (out/'STAR_Log.out').read_text()
   qc={k.strip():v.strip() for l in (out/'STAR_Log.final.out').read_text().splitlines() if '|' in l for k,v in [l.split('|',1)]}
   data=[l.split('\t') for l in (out/'STAR_ReadsPerGene.out.tab').read_text().splitlines()];assert int(qc['Number of input reads'])==N and [v[0] for v in data[4:]]==expected
   assert all(sum(int(v[i]) for v in data)==N for i in [1,2,3]);vec=np.array([int(v[1]) for v in data[4:]]);vectors[mode]=vec
   h=sha(out/'STAR_ReadsPerGene.out.tab')
   if receipt.exists():assert json.loads(receipt.read_text())['counts_sha256']==h
   else:write(receipt,dict(signature=signature,counts_sha256=h,QC=qc,completed=time.time()))
   metrics.append(dict(sample=r['column'],run=run,mode=mode,input_units=N,unique_pct=float(qc['Uniquely mapped reads %'].strip('%')),assigned_unstranded=int(vec.sum()),assigned_pct=float(vec.sum()*100/N),detected_genes=int((vec>0).sum())))
   state['completed'].append(run+'/'+mode);save()
  for a,b in [('R1','paired'),('R2','paired'),('R1','R2')]:
   mask=(vectors[a]+vectors[b])>0
   comparisons.append(dict(sample=r['column'],mode_a=a,mode_b=b,genes_nonzero_either=int(mask.sum()),spearman_nonzero_either=float(spearmanr(vectors[a][mask],vectors[b][mask]).statistic)))
 for filename,data in [('mate_alignment_QC.tsv',metrics),('mate_count_concordance.tsv',comparisons)]:
  with (O/filename).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
 write(O/'execution_complete.json',dict(samples=9,alignments=27,subset_pairs_per_sample=N,metrics=metrics,comparisons=comparisons,limitations=['Initial subset, not randomized or full-library verification','Single-end input units are reads; paired input units are pairs','No proof of library kit, physical sample identity or biological safety','Different sensitivity/multimapping of single versus paired alignment expected','Subset too shallow to interpret zero counts of low-expression candidate-neighbor genes']))
 state['status']='completed_requires_review';save()
except Exception as e:state.update(status='failed',error=str(e),traceback=traceback.format_exc());save();raise
