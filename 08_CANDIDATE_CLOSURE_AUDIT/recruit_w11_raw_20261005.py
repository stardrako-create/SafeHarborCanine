from pathlib import Path
import pysam,subprocess,json,datetime,fcntl,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_RAW_RECRUIT_20261005';O.mkdir(exist_ok=False);lock=open(O/'worker.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def state(stage,**kw):(O/'status.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':stage,**kw},indent=2))
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with pysam.FastaFile(str(A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa')) as f:j=f.fetch('junction_L4707939_R33613308').upper()
with pysam.FastaFile(str(A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa')) as f:other=rc(f.fetch('NC_051836.1',35342033,35342133).upper())
seeds={'junction':j[980:1020],'tail':'A'*20+other[:20],'left_flank':j[950:1000],'opposite_flank':other[:50]};(O/'seeds.json').write_text(json.dumps(seeds,indent=2));patterns=[(k,s.encode()) for k,v in seeds.items() for s in [v,rc(v)]]
cmd=['/usr/bin/fastq-dump-orig','--split-spot','--stdout','-A','SRR14750332','/mnt/c/Users/Utilizador/SafeHarborData/SRR14750332/SRR14750332.sra'];(O/'command.json').write_text(json.dumps(cmd));n=kept=0
try:
 state('scanning_raw_spots',pairs=0)
 with (O/'extract.log').open('wb') as err,(O/'selected_R1.fastq').open('wb') as a,(O/'selected_R2.fastq').open('wb') as b,(O/'matches.jsonl').open('w') as matches:
  proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err)
  while True:
   x=[proc.stdout.readline() for _ in range(4)]
   if not x[0]:break
   y=[proc.stdout.readline() for _ in range(4)]
   assert x[0].startswith(b'@') and x[0].split()[0]==y[0].split()[0]
   for r in [x,y]:assert r[2].startswith(b'+') and len(r[1].strip())==len(r[3].strip())==150
   n+=1;found=sorted({k for k,p in patterns if p in x[1] or p in y[1]})
   if found:
    kept+=1;a.writelines(x);b.writelines(y);matches.write(json.dumps({'name':x[0].split()[0][1:].decode(),'seeds':found})+'\n')
   if n%1000000==0:state('scanning_raw_spots',pairs=n,selected_pairs=kept)
  assert proc.wait()==0 and n==214126885
 state('completed_needs_review',pairs=n,selected_pairs=kept,sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [O/'selected_R1.fastq',O/'selected_R2.fastq',O/'matches.jsonl']})
except Exception as e:
 if 'proc' in globals() and proc.poll() is None:proc.terminate()
 state('failed',pairs=n,selected_pairs=kept,error=str(e));raise