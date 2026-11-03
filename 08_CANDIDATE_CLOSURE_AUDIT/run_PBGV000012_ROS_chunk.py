import pathlib,subprocess,json,shutil,fcntl,datetime,sys
BASE=pathlib.Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
OUT=BASE/'08_CANDIDATE_CLOSURE_AUDIT/SV_PBGV000012/ROS_chunked';OUT.mkdir(parents=True,exist_ok=True)
lock=open(OUT/'worker.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
BIN='/home/stardrako/miniforge3/envs/atac/bin/'
SRA='/mnt/c/Users/Utilizador/SafeHarborData/SRR14750332/SRR14750332.sra'
REF=BASE/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'
start,end=1,1000000
D=OUT/'chunk_000001_1000000';D.mkdir(exist_ok=True)
assert not (D/'receipt.json').exists()
def run(args,name):
 with open(D/(name+'.out'),'w') as o,open(D/(name+'.err'),'w') as e:subprocess.run(list(map(str,args)),stdout=o,stderr=e,check=True)
def state(x): (OUT/'status.json').write_text(json.dumps({'stage':x,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'start':start,'end':end}))
try:
 assert shutil.disk_usage(OUT).free>5*1024**3
 state('extracting')
 run(['/usr/bin/fastq-dump-orig','--split-files','-N',start,'-X',end,'-A','SRR14750332','-O',D,SRA],'extract')
 fq=[D/f'SRR14750332_{i}.fastq' for i in (1,2)]
 with open(fq[0]) as f,open(fq[1]) as g:
  n=0
  while True:
   a=[f.readline().rstrip() for _ in range(4)];b=[g.readline().rstrip() for _ in range(4)]
   if not a[0] and not b[0]:break
   assert a[0].split()[0]==b[0].split()[0]
   for x in (a,b):assert x[0].startswith('@') and x[2].startswith('+') and len(x[1])==len(x[3])==150
   n+=1
  assert n==end-start+1
 state('aligning_full_ROS_retaining_chr32_links')
 import pysam
 with open(D/'bwa.err','w') as e:
  proc=subprocess.Popen([BIN+'bwa','mem','-t','4','-R','@RG\tID:SRR14750332\tSM:PBGV000012\tPL:ILLUMINA',str(REF),*map(str,fq)],stdout=subprocess.PIPE,stderr=e)
  total=primary=kept=0
  with pysam.AlignmentFile(proc.stdout,'r') as sam:
   tid=sam.get_tid('NC_051836.1');assert tid>=0
   with pysam.AlignmentFile(str(D/'selected.unsorted.bam'),'wb',template=sam) as bam:
    for r in sam:
     total+=1;primary+=not(r.is_secondary or r.is_supplementary)
     if r.reference_id==tid or r.next_reference_id==tid:bam.write(r);kept+=1
  assert proc.wait()==0 and primary==2*n
 run([BIN+'samtools','sort','-@','2','-m','256M','-o',D/'selected.bam',D/'selected.unsorted.bam'],'sort')
 run([BIN+'samtools','quickcheck','-v',D/'selected.bam'],'quickcheck')
 run([BIN+'samtools','index',D/'selected.bam'],'index')
 (D/'receipt.json').write_text(json.dumps({'start':start,'end':end,'pairs':n,'alignment_records':total,'primary_records':primary,'selected_records':kept,'reference':str(REF),'selection':'reference or mate on NC_051836.1','duplicates':'not yet marked','scope':'first million spots only; no structural conclusion'},indent=2))
 state('pilot_chunk_complete')
except Exception as e:
 state('failed: '+str(e));raise