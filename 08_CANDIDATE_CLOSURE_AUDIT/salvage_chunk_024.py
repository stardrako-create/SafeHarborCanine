from pathlib import Path
import sys, os, json, time, fcntl, subprocess
import pysam
from recover_sv_chunks import A, B, CHROM, sha, atomic
lock=(A/'worker.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
c=A/'chunk_024'; receipt=c/'validated.json'; assert not receipt.exists()
state=A/'salvage_chunk_024_status.json'
def status(stage, **kw):
 atomic(state,dict(stage=stage,pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**kw))
try:
 status('validating_complete_SAM')
 assert '[main] Real time:' in (A/'chunk_024_align.stderr.log').read_text()
 dest=c/'retained_recovered.bam'; assert not dest.exists()
 primary=r1=r2=records=retained=0
 with pysam.AlignmentFile(str(c/'replay.sam'),'r') as src, pysam.AlignmentFile(str(dest),'wb',template=src) as out:
  tid=src.get_tid(CHROM);assert tid>=0
  for r in src:
   records+=1
   if not r.is_secondary and not r.is_supplementary:
    primary+=1;r1+=int(r.is_read1);r2+=int(r.is_read2)
   if r.reference_id==tid or r.next_reference_id==tid:out.write(r);retained+=1
   if records%1000000==0:status('validating_complete_SAM',records=records)
 assert (primary,r1,r2)==(8000000,4000000,4000000),(primary,r1,r2)
 subprocess.run([str(B/'samtools'),'quickcheck','-v',str(dest)],check=True)
 old=c/'retained.bam'; saved=c/'retained.interrupted_20260921T1513.bam'; assert not saved.exists()
 if old.exists():old.rename(saved)
 dest.rename(old)
 atomic(receipt,dict(start_pair0=96000000,pairs=4000000,primary=primary,read1=r1,read2=r2,all_records=records,retained_records=retained,bam_sha256=sha(old),validated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),recovery='Completed SAM re-read to EOF; interrupted BAM preserved separately'))
 status('validated_resuming_main_worker',pairs=4000000)
except BaseException as e:
 status('failed_preserved',error=repr(e));raise
finally:
 fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
os.execv(sys.executable,[sys.executable,str(Path(__file__).with_name('recover_sv_chunks.py'))])
