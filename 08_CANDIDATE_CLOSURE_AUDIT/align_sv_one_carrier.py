"""Full-genome mapping, retained chr32/mate alignments, disk-guarded diagnostics."""
from pathlib import Path
import subprocess,json,os,time,fcntl,shutil,signal
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT';D=O/'SV_PBGV000010';A=D/'alignment';BIN=Path('/home/stardrako/miniforge3/envs/atac/bin')
REF=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');CHROM='NC_049253.1'

def main():
 qc=json.loads((D/'fastq_pair_validation.json').read_text());assert qc['status']=='validated_FASTQ_structure_and_mate_identity' and qc['gzip_streams_read_to_EOF']
 A.mkdir(exist_ok=True);lock=(A/'worker.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 assert not (A/'completed.json').exists(),'Completed outputs already exist'
 assert shutil.disk_usage(A).free>30*1024**3,'Need at least30GiB free before alignment'
 # Validate receipts and reference metadata without changing input files.
 for i in [1,2]:
  p=D/f'SRR15734832_{i}.fastq.gz';r=json.loads((D/(p.name+'.validated.json')).read_text());assert r['validated'] and p.stat().st_size==r['bytes']
 refcheck=json.loads((O/'UU_reference_preflight.json').read_text());assert refcheck['BWA_ann_matches_FAI_all_contigs']
 state={'pid':os.getpid(),'status':'starting','stage':'initialization','input_pairs':qc['pairs']}
 def save(**kw):
  state.update(kw);state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());state['free_bytes']=shutil.disk_usage(A).free
  t=A/'status.tmp';t.write_text(json.dumps(state,indent=2));os.replace(t,A/'status.json')
 processes=[];logs=[]
 def monitor(ps):
  while any(p.poll() is None for p in ps):
   if shutil.disk_usage(A).free<20*1024**3:raise RuntimeError('20GiB disk floor reached; preserve all partial outputs')
   if any(p.poll() not in [None,0] for p in ps):raise RuntimeError('Pipeline process failed; inspect stage stderr')
   save(child_pids=[p.pid for p in ps]);time.sleep(5)
  assert all(p.returncode==0 for p in ps),'Nonzero pipeline exit'
 def run(cmd,label):
  save(stage=label,status='running');log=(A/(label+'.stderr.log')).open('ab');logs.append(log)
  p=subprocess.Popen(cmd,stderr=log,stdout=subprocess.DEVNULL,start_new_session=True);processes.append(p);monitor([p])
 try:
  # Full reference participates in alignment competition; only storage is regional.
  commands=[
   [str(BIN/'bwa'),'mem','-t','8','-R','@RG\\tID:PBGV000010\\tSM:PBGV000010\\tPL:ILLUMINA',str(REF),str(D/'SRR15734832_1.fastq.gz'),str(D/'SRR15734832_2.fastq.gz')],
   [str(BIN/'samtools'),'view','-u','--save-counts',str(A/'selection_counts.json'),'-e',f'rname == "{CHROM}" || rnext == "{CHROM}"','-'],
   [str(BIN/'samtools'),'sort','-n','-@','2','-m','512M','-T',str(A/'name_sort_tmp'),'-o',str(A/'selected_namesort.bam'),'-']]
  (A/'commands.json').write_text(json.dumps(commands,indent=2));save(stage='full_genome_alignment_regional_retention',status='running')
  previous=None;ps=[]
  for i,cmd in enumerate(commands):
   log=(A/f'pipeline_{i}.stderr.log').open('ab');logs.append(log)
   p=subprocess.Popen(cmd,stdin=previous,stdout=subprocess.PIPE if i<len(commands)-1 else subprocess.DEVNULL,stderr=log,start_new_session=True)
   if previous is not None:previous.close()
   previous=p.stdout;ps.append(p);processes.append(p)
  monitor(ps)
  run([str(BIN/'samtools'),'fixmate','-m',str(A/'selected_namesort.bam'),str(A/'selected_fixmate.bam')],'fixmate')
  run([str(BIN/'samtools'),'sort','-@','2','-m','512M','-T',str(A/'coord_sort_tmp'),'-o',str(A/'selected_coordsort.bam'),str(A/'selected_fixmate.bam')],'coordinate_sort')
  run([str(BIN/'samtools'),'markdup','-s','-f',str(A/'markdup_stats.txt'),str(A/'selected_coordsort.bam'),str(A/'selected_markdup.bam')],'mark_duplicates')
  run([str(BIN/'samtools'),'quickcheck','-v',str(A/'selected_markdup.bam')],'quickcheck')
  run([str(BIN/'samtools'),'index',str(A/'selected_markdup.bam')],'index')
  # Retain all duplicates; analysis can compare inclusion/exclusion using the flag.
  out={'status':'alignment_completed_pending_biological_review','input_pairs':qc['pairs'],'reference':str(REF),'retention':CHROM+' or mate on this chromosome; full genome used for mapping','bam':str(A/'selected_markdup.bam'),'limitations':['Regional retention omits unrelated genome alignments','Duplicate marking on retained subset is not equivalent to whole-genome deduplication for all cross-chromosome fragments','No sample identity or SV truth established by successful alignment','All intermediates retained; no cleanup performed']}
  (A/'completed.json').write_text(json.dumps(out,indent=2));save(status='completed',stage='awaiting_coverage_and_breakpoint_review')
 except BaseException as e:
  for p in processes:
   if p.poll() is None:
    try:os.killpg(p.pid,signal.SIGTERM)
    except ProcessLookupError:pass
  for p in processes:
   try:p.wait(timeout=15)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
  save(status='failed',error=repr(e));raise
 finally:
  for log in logs:log.close()

if __name__=='__main__':main()
