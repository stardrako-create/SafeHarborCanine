"""Restartable full-reference WGS alignment; bounded temporary files and receipts."""
from pathlib import Path
import json, os, time, gzip, shutil, subprocess, signal, fcntl, hashlib
import pysam

P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT'; D=O/'SV_PBGV000010'; A=D/'alignment_chunked'
B=Path('/home/stardrako/miniforge3/envs/atac/bin'); REF=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa')
CHROM='NC_049253.1'; SIZE=4_000_000

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def atomic(p,value):
 t=p.with_suffix('.tmp');t.write_text(json.dumps(value,indent=2));os.replace(t,p)

def main():
 A.mkdir(exist_ok=True);lock=(A/'worker.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 state={'pid':os.getpid()}
 def status(stage,**kw):
  state.update(stage=stage,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**kw)
  atomic(A/'status.json',state)
 def disk():
  if shutil.disk_usage(A).free<20*1024**3:raise RuntimeError('Below 20 GiB disk floor; partial files preserved')
 def run(cmd,label,out=None):
  disk();status(label)
  with (A/(label+'.stderr.log')).open('ab') as err:
   stdout=out.open('wb') if out else open(os.devnull,'wb')
   with stdout:
    proc=subprocess.Popen(cmd,stdout=stdout,stderr=err,start_new_session=True)
    try:
     while proc.poll() is None:
      disk();status(label,child_pid=proc.pid);time.sleep(5)
     if proc.returncode:raise RuntimeError(f'{label}: exit {proc.returncode}')
    except BaseException:
     if proc.poll() is None:
      os.killpg(proc.pid,signal.SIGTERM)
      try:proc.wait(timeout=15)
      except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);proc.wait()
     raise
 def clean_generated(c):
  assert c.resolve().parent==A.resolve()
  for name in ('mate1.fastq','mate2.fastq','replay.sam'):
   f=c/name
   assert f.resolve().parent==c.resolve()
   if f.exists():f.unlink()
 try:
  qc=json.loads((D/'fastq_pair_validation.json').read_text());total=qc['pairs']
  assert qc['gzip_streams_read_to_EOF'] and qc['status']=='validated_FASTQ_structure_and_mate_identity'
  for mate in (1,2):
   f=D/f'SRR15734832_{mate}.fastq.gz';r=json.loads((D/(f.name+'.validated.json')).read_text())
   assert r['validated'] and f.stat().st_size==r['bytes']
  assert not (A/'completed.json').exists()
  assert shutil.disk_usage(A).free>30*1024**3
  atomic(A/'method.json',{'chunk_pairs':SIZE,'input_pairs':total,'reference':str(REF),'retention':CHROM+' or mate on chromosome','notes':['Full reference used in each chunk; insert-size estimation is per chunk, not identical to original monolithic run','All chunks combined before fixmate and duplicate marking','Only generated FASTQ/SAM of receipt-validated chunks removed; original inputs and failed run preserved']})
  bams=[]
  with gzip.open(D/'SRR15734832_1.fastq.gz','rb') as f1,gzip.open(D/'SRR15734832_2.fastq.gz','rb') as f2:
   for index,start in enumerate(range(0,total,SIZE)):
    n=min(SIZE,total-start);c=A/f'chunk_{index:03d}';c.mkdir(exist_ok=True);receipt=c/'validated.json';bam=c/'retained.bam'
    valid=receipt.exists()
    if valid:
     r=json.loads(receipt.read_text());assert r['start_pair0']==start and r['pairs']==n and sha(bam)==r['bam_sha256']
     clean_generated(c)
    else:
     assert not any((c/x).exists() for x in ('mate1.fastq','mate2.fastq','replay.sam','retained.bam')),'Unvalidated partial chunk exists; preserve and investigate before resuming'
    status('extracting_chunk' if not valid else 'skipping_validated_chunk',chunk=index,start_pair0=start,pairs=n)
    handles=[] if valid else [(c/f'mate{i}.fastq').open('wb') for i in (1,2)]
    try:
     for j in range(n):
      r1=[f1.readline() for _ in range(4)];r2=[f2.readline() for _ in range(4)]
      assert r1[0].startswith(b'@') and r2[0].startswith(b'@') and r1[2].startswith(b'+') and r2[2].startswith(b'+') and r1[3] and r2[3]
      assert r1[0].split()[0].removesuffix(b'/1')==r2[0].split()[0].removesuffix(b'/2')
      if not valid:handles[0].writelines(r1);handles[1].writelines(r2)
      if j%100_000==0:disk();status('skipping_validated_chunk' if valid else 'extracting_chunk',chunk_pairs_read=j)
    finally:
     for h in handles:h.close()
    if not valid:
     cmd=[str(B/'bwa'),'mem','-t','8','-R','@RG\\tID:PBGV000010\\tSM:PBGV000010\\tPL:ILLUMINA',str(REF),str(c/'mate1.fastq'),str(c/'mate2.fastq')]
     atomic(c/'command.json',cmd);run(cmd,f'chunk_{index:03d}_align',c/'replay.sam')
     primary=read1=read2=records=retained=0
     status('validate_and_retain',chunk=index)
     with pysam.AlignmentFile(str(c/'replay.sam'),'r') as src,pysam.AlignmentFile(str(bam),'wb',template=src) as out:
      tid=src.get_tid(CHROM);assert tid>=0
      for r in src:
       records+=1
       if not r.is_secondary and not r.is_supplementary:
        primary+=1;read1+=int(r.is_read1);read2+=int(r.is_read2)
       if r.reference_id==tid or r.next_reference_id==tid:out.write(r);retained+=1
       if records%100_000==0:disk();status('validate_and_retain',records_read=records)
     assert primary==2*n and read1==n and read2==n,(primary,read1,read2,n)
     run([str(B/'samtools'),'quickcheck','-v',str(bam)],f'chunk_{index:03d}_quickcheck')
     atomic(receipt,{'start_pair0':start,'pairs':n,'primary':primary,'read1':read1,'read2':read2,'all_records':records,'retained_records':retained,'bam_sha256':sha(bam),'validated_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})
     clean_generated(c)
    bams.append(bam)
   assert not f1.read(1) and not f2.read(1),'Unexpected extra FASTQ data'
  atomic(A/'all_chunks_validated.json',{'pairs':total,'chunks':len(bams)})
  listing=A/'bam_list.txt';listing.write_text('\n'.join(map(str,bams))+'\n')
  run([str(B/'samtools'),'cat','-b',str(listing),'-o',str(A/'retained_all.bam')],'concatenate')
  run([str(B/'samtools'),'sort','-n','-@','2','-m','512M','-o',str(A/'selected_namesort.bam'),str(A/'retained_all.bam')],'namesort')
  run([str(B/'samtools'),'fixmate','-m',str(A/'selected_namesort.bam'),str(A/'selected_fixmate.bam')],'fixmate')
  run([str(B/'samtools'),'sort','-@','2','-m','512M','-o',str(A/'selected_coordsort.bam'),str(A/'selected_fixmate.bam')],'coordsort')
  run([str(B/'samtools'),'markdup','-s','-f',str(A/'markdup_stats.txt'),str(A/'selected_coordsort.bam'),str(A/'selected_markdup.bam')],'markdup')
  run([str(B/'samtools'),'quickcheck','-v',str(A/'selected_markdup.bam')],'final_quickcheck')
  run([str(B/'samtools'),'index',str(A/'selected_markdup.bam')],'index')
  atomic(A/'completed.json',{'status':'alignment_completed_pending_biological_review','input_pairs':total,'chunks':len(bams),'bam':str(A/'selected_markdup.bam'),'limitations':['Regional duplicate marking','Chunk-specific insert size estimation','No biological or SV validation implied']})
  status('completed_pending_review',child_pid=None)
 except BaseException as e:status('failed',error=repr(e));raise

if __name__=='__main__':main()
