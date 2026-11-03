import pathlib,subprocess,json,datetime,fcntl,hashlib,shutil
A=pathlib.Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/08_CANDIDATE_CLOSURE_AUDIT')
O=A/'SV_PBGV000012/ROS_production_20260930';O.mkdir(exist_ok=True)
l=open(O/'supervisor.lock','w');fcntl.flock(l,fcntl.LOCK_EX|fcntl.LOCK_NB)
TOTAL=214126885
PY='/home/stardrako/miniforge3/envs/atac/bin/python'
def state(s,**kw):
 (O/'supervisor_status.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':s,'total_spots':TOTAL,**kw},indent=2))
try:
 first=A/'SV_PBGV000012/ROS_chunked_v2/chunk_000001_1000000/receipt.json'
 assert json.loads(first.read_text())['primary_records']==2000000
 for start in range(1000001,TOTAL+1,1000000):
  end=min(start+999999,TOTAL);d=O/f'chunk_{start:09d}_{end:09d}'
  if (d/'receipt.json').exists():continue
  assert shutil.disk_usage(O).free>5*1024**3
  state('processing',start=start,end=end)
  with open(O/f'launcher_{start}.log','a') as log:
   subprocess.run([PY,str(A/'run_PBGV000012_ROS_production_chunk.py'),str(start),str(end)],stdout=log,stderr=log,check=True)
  r=json.loads((d/'receipt.json').read_text());assert r['pairs']==end-start+1 and r['primary_records']==2*r['pairs']
  hashes={}
  for f in [d/'selected.bam',d/'selected.bam.bai',d/'SRR14750332_1.fastq',d/'SRR14750332_2.fastq']:
   h=hashlib.sha256()
   with f.open('rb') as stream:
    for b in iter(lambda:stream.read(8*1024**2),b''):h.update(b)
   hashes[f.name]=h.hexdigest()
  (d/'sha256.json').write_text(json.dumps(hashes,indent=2))
  # Only reproducible extraction intermediates; retain original SRA, all BAMs and receipts.
  for name in ['SRR14750332_1.fastq','SRR14750332_2.fastq']:
   f=d/name;assert f.resolve().parent==d.resolve() and not f.is_symlink();f.unlink()
  state('chunk_validated',completed_through=end)
 state('all_chunks_completed_not_merged_or_deduplicated',completed_through=TOTAL)
except Exception as e:state('failed',error=str(e));raise