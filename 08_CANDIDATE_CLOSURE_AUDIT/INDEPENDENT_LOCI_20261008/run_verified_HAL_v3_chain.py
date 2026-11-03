from pathlib import Path
import json,subprocess,datetime,os,fcntl,hashlib,traceback
N=Path(__file__).resolve().parent;O=N/'verified_HAL_v3_chain_20261009';O.mkdir(exist_ok=False);lock=open('/tmp/safeharbor_hal_verified_chain.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);python='/home/stardrako/miniforge3/envs/atac/bin/python'
status=dict(state='running',pid=os.getpid(),utc_start=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='extracting',scripts_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [N/'extract_verified_HAL_v3.py',N/'score_verified_MAF_v3.py']})
def save():
 tmp=O/'status.tmp'
 with tmp.open('w') as f:json.dump(status,f,indent=2);f.flush();os.fsync(f.fileno())
 tmp.replace(O/'status.json')
try:
 save()
 with (O/'extraction.out').open('w') as out,(O/'extraction.err').open('w') as err:subprocess.run([python,str(N/'extract_verified_HAL_v3.py')],stdout=out,stderr=err,check=True)
 receipt=json.loads((N/'remote_HAL_verified_20261009_v3/status.json').read_text());assert receipt['state']=='completed','Extraction not validated; scoring forbidden';status['stage']='scoring';save()
 with (O/'scoring.out').open('w') as out,(O/'scoring.err').open('w') as err:subprocess.run([python,str(N/'score_verified_MAF_v3.py')],stdout=out,stderr=err,check=True)
 receipt=json.loads((N/'remote_conservation_verified_20261009_v3/status.json').read_text());assert receipt['state']=='completed' and len(receipt['results'])==3 and not receipt['errors'];status.update(state='completed',stage='validated',results=receipt['results'])
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:status['utc_end']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
