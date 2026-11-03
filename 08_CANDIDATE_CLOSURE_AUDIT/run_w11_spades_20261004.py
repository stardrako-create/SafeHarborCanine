from pathlib import Path
import tarfile,subprocess,json,hashlib,datetime,fcntl
A=Path(__file__).resolve().parent;O=A/'W11_SPADES_20261004';lock=open(O/'worker.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def state(s,**kw):(O/'status.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':s,**kw},indent=2))
try:
 archive=O/'SPAdes-4.3.0-Linux.tar.gz';state('extracting',archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest())
 with tarfile.open(archive) as t:
  for m in t.getmembers():
   assert (O/m.name).resolve().is_relative_to(O.resolve())
  t.extractall(O,filter='data')
 exe=O/'SPAdes-4.3.0-Linux/bin/spades.py'
 version=subprocess.check_output(['python3',str(exe),'--version'],text=True);(O/'version.txt').write_text(version)
 cmd=['python3',str(exe),'--only-assembler','-1',str(O/'R1.fastq'),'-2',str(O/'R2.fastq'),'-s',str(O/'single.fastq'),'-k','21,33,55,77','-t','6','-m','12','-o',str(O/'assembly')]
 (O/'command.json').write_text(json.dumps(cmd));state('assembling',version=version.strip())
 with (O/'run.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 assert (O/'assembly/contigs.fasta').exists();state('assembly_completed_needs_review',contigs_sha256=hashlib.sha256((O/'assembly/contigs.fasta').read_bytes()).hexdigest())
except Exception as e:state('failed',error=str(e));raise