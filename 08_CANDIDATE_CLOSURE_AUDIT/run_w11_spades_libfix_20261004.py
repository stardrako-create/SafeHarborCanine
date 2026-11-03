from pathlib import Path
import os,json,subprocess,datetime,hashlib,fcntl
A=Path(__file__).resolve().parent;O=A/'W11_SPADES_20261004';W=Path('/mnt/d/SafeHarborWork/w11_spades_20261004_ascii');lock=open(W/'worker.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def state(s,**kw):(O/'status_libfix.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':s,**kw},indent=2))
try:
 env=os.environ.copy();env['LD_LIBRARY_PATH']='/home/stardrako/miniforge3/envs/atac/lib'
 deps=subprocess.check_output(['ldd',str(W/'SPAdes-4.3.0-Linux/bin/spades-core')],env=env,text=True);(O/'ldd_libfix.txt').write_text(deps);assert 'not found' not in deps
 cmd=json.loads((O/'command_ascii.json').read_text());cmd[-1]=str(W/'assembly_libfix');assert not Path(cmd[-1]).exists();(O/'command_libfix.json').write_text(json.dumps(cmd));state('assembling',LD_LIBRARY_PATH=env['LD_LIBRARY_PATH'])
 with (W/'run_libfix.log').open('w') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 f=W/'assembly_libfix/contigs.fasta';assert f.exists();state('assembly_completed_needs_review',contigs=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest())
except Exception as e:state('failed',error=str(e));raise