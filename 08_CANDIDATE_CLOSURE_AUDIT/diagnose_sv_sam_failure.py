"""Replay a bounded input neighborhood; preserve SAM to inspect parser failures."""
from pathlib import Path
import gzip,json,time,subprocess,fcntl,shutil,os,signal,re
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT';D=O/'SV_PBGV000010';T=D/'sam_failure_diagnostic';B=Path('/home/stardrako/miniforge3/envs/atac/bin')
START=106_000_000;STOP=110_000_000
def main():
 T.mkdir(exist_ok=True);lock=(T/'worker.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 def status(stage,**kw):
  t=T/'status.tmp';t.write_text(json.dumps(dict(stage=stage,pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**kw),indent=2));os.replace(t,T/'status.json')
 try:
  assert shutil.disk_usage(T).free>25*1024**3
  for mate in [1,2]:
   dest=T/f'mate{mate}.fastq.gz'
   assert not dest.exists(),'Existing diagnostic input preserved; do not overwrite'
   with gzip.open(D/f'SRR15734832_{mate}.fastq.gz','rb') as src,gzip.open(dest,'wb',compresslevel=1) as out:
    for i in range(STOP):
     rec=[src.readline() for _ in range(4)]
     assert rec[0].startswith(b'@') and rec[2].startswith(b'+') and rec[3],f'Truncated input at {i}'
     if i>=START:out.writelines(rec)
     if i%1000000==0:status('extracting_bounded_replay',mate=mate,scanned_pairs=i,first_pair0=START,last_pair0_exclusive=STOP)
     if i%1000000==0 and shutil.disk_usage(T).free<20*1024**3:raise RuntimeError('Disk floor reached')
  def run(cmd,stage,out):
   status(stage,command=cmd)
   with (T/(stage+'.stderr.log')).open('wb') as err,out.open('wb') as stdout:
    p=subprocess.Popen(cmd,stdout=stdout,stderr=err,start_new_session=True)
    try:
     while p.poll() is None:
      if shutil.disk_usage(T).free<20*1024**3:raise RuntimeError('Disk floor reached')
      status(stage,child_pid=p.pid,output_bytes=out.stat().st_size);time.sleep(5)
    except BaseException:
     os.killpg(p.pid,signal.SIGTERM);p.wait();raise
    return p.returncode
  cmd=[str(B/'bwa'),'mem','-t','8','-R','@RG\\tID:PBGV000010\\tSM:PBGV000010\\tPL:ILLUMINA','/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa',str(T/'mate1.fastq.gz'),str(T/'mate2.fastq.gz')]
  rc=run(cmd,'bwa_replay',T/'replay.sam');assert rc==0,'BWA replay failed'
  rc=run([str(B/'samtools'),'view','-b','-o',str(T/'replay.bam'),str(T/'replay.sam')],'sam_parse',T/'sam_parse.stdout.log')
  if rc:
   text=(T/'sam_parse.stderr.log').read_text();matches=re.findall(r'line (\d+)',text)
   if matches:
    target=int(matches[-1])
    with (T/'replay.sam').open('rb') as src,(T/'failing_line_context.sam.txt').open('wb') as out:
     for i,line in enumerate(src,1):
      if target-2<=i<=target+2:out.write(str(i).encode()+b'\t'+line)
      if i>target+2:break
  status('replay_complete_pending_review',sam_parse_exit_code=rc,scope='Bounded diagnostic; original failed full alignment remains invalid; successful replay does not establish original failure cause')
 except BaseException as e:status('failed',error=repr(e));raise
if __name__=='__main__':main()
