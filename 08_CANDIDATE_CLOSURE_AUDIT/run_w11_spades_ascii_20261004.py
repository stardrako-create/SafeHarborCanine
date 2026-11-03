from pathlib import Path
import tarfile,subprocess,json,hashlib,datetime,shutil,fcntl
A=Path(__file__).resolve().parent;O=A/'W11_SPADES_20261004';W=Path('/mnt/d/SafeHarborWork/w11_spades_20261004_ascii');W.mkdir(parents=True,exist_ok=False)
def state(s,**kw):(O/'status_ascii.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':s,'workspace':str(W),**kw},indent=2))
lock=open(W/'worker.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
 state('preparing_ASCII_workspace')
 for name in ['R1.fastq','R2.fastq','single.fastq']:
  shutil.copy2(O/name,W/name);assert hashlib.sha256((O/name).read_bytes()).digest()==hashlib.sha256((W/name).read_bytes()).digest()
 with tarfile.open(O/'SPAdes-4.3.0-Linux.tar.gz') as t:t.extractall(W,filter='data')
 exe=W/'SPAdes-4.3.0-Linux/bin/spades.py';cmd=['python3',str(exe),'--only-assembler','-1',str(W/'R1.fastq'),'-2',str(W/'R2.fastq'),'-s',str(W/'single.fastq'),'-k','21,33,55,77','-t','6','-m','12','-o',str(W/'assembly')]
 (O/'command_ascii.json').write_text(json.dumps(cmd));state('assembling')
 with (W/'run.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 contigs=W/'assembly/contigs.fasta';assert contigs.exists();shutil.copy2(contigs,O/'contigs_ascii.fasta');state('assembly_completed_needs_review',contigs_sha256=hashlib.sha256(contigs.read_bytes()).hexdigest())
except Exception as e:state('failed',error=str(e));raise