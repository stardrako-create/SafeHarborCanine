from pathlib import Path
import pysam,json,datetime,subprocess,os,hashlib,fcntl
A=Path(__file__).resolve().parent;O=A/'W11_SPADES_PBGV10_20261005';O.mkdir(exist_ok=False);W=Path('/mnt/d/SafeHarborWork/w11_PBGV10_20261005');W.mkdir(exist_ok=False);l=open(W/'worker.lock','w');fcntl.flock(l,fcntl.LOCK_EX|fcntl.LOCK_NB)
def state(stage,**kw):(O/'status.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':stage,**kw},indent=2))
try:
 state('extracting_local_reads');B=A/'SV_PBGV000010/alignment_chunked/selected_markdup.bam';names=set();regions=[(4687939,4727939),(33583308,33643308)]
 with pysam.AlignmentFile(str(B),'rb') as f:
  for s,e in regions:
   for r in f.fetch('NC_049253.1',s,e):
    if not(r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail):names.add(r.query_name)
 reads={}
 with pysam.AlignmentFile(str(B),'rb') as f:
  for r in f:
   if r.query_name in names and not(r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail):
    m=1 if r.is_read1 else 2;assert m not in reads.get(r.query_name,{});reads.setdefault(r.query_name,{})[m]=(r.get_forward_sequence(),r.get_forward_qualities())
 n=single=0
 with (W/'R1.fastq').open('w') as a,(W/'R2.fastq').open('w') as b,(W/'single.fastq').open('w') as c:
  for name,rr in sorted(reads.items()):
   if len(rr)==2:n+=1
   else:single+=1
   for m,(seq,q) in rr.items():
    f=([a,b][m-1] if len(rr)==2 else c);assert len(seq)==len(q);f.write('@'+name+'/'+str(m)+'\n'+seq+'\n+\n'+''.join(chr(x+33) for x in q)+'\n')
 (O/'inputs.json').write_text(json.dumps({'pairs':n,'singles':single,'regions_UU0':regions,'selection':'primary nonduplicate QCpass, no MAPQ threshold, mates retained in BAM','sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in W.glob('*.fastq')}},indent=2))
 exe='/mnt/d/SafeHarborWork/w11_spades_20261004_ascii/SPAdes-4.3.0-Linux/bin/spades.py';cmd=['python3',exe,'--only-assembler','-1',str(W/'R1.fastq'),'-2',str(W/'R2.fastq'),'-k','21,33,55,77','-t','6','-m','12','-o',str(W/'assembly')]
 if single:cmd+=['-s',str(W/'single.fastq')]
 env=os.environ.copy();env['LD_LIBRARY_PATH']='/home/stardrako/miniforge3/envs/atac/lib';(O/'command.json').write_text(json.dumps(cmd));state('assembling',pairs=n,singles=single)
 with (W/'run.log').open('w') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 c=W/'assembly/contigs.fasta';assert c.exists();state('completed_needs_review',contigs=str(c),sha256=hashlib.sha256(c.read_bytes()).hexdigest())
except Exception as e:state('failed',error=str(e));raise