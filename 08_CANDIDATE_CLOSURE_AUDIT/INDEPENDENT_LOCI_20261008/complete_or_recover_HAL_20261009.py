from pathlib import Path
import os,select,json,datetime,subprocess,traceback,fcntl,hashlib
N=Path(__file__).resolve().parent;O=N/'HAL_completion_and_recovery_20261009';O.mkdir(exist_ok=False)
lock=open('/tmp/safeharbor_HAL_completion_recovery.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
status=dict(state='running',pid=os.getpid(),stage='waiting_for_existing_v3_chain',utc_start=datetime.datetime.now(datetime.timezone.utc).isoformat(),watched_pid=292)
def save():
 tmp=O/'status.tmp';tmp.write_text(json.dumps(status,indent=2));tmp.replace(O/'status.json')
python='/home/stardrako/miniforge3/envs/atac/bin/python'
def run(script,label):
 status.update(stage=label,script=script,script_sha256=hashlib.sha256((N/script).read_bytes()).hexdigest());save()
 with (O/(label+'.out')).open('w') as out,(O/(label+'.err')).open('w') as err:subprocess.run([python,str(N/script)],stdout=out,stderr=err,check=True)
try:
 save()
 # OS wait on existing process termination: no periodic polling or duplicate extraction.
 try:
  fd=os.pidfd_open(292);select.select([fd],[],[]);os.close(fd)
 except ProcessLookupError:pass
 old=json.loads((N/'verified_HAL_v3_chain_20261009/status.json').read_text());status['v3_terminal_receipt']=old
 if old['state']=='completed':run('review_remote_conservation_taxa.py','review_v3_taxa');status.update(state='completed',version='v3',stage='validated_taxa_review')
 elif old['state']=='failed':
  extract=json.loads((N/'remote_HAL_verified_20261009_v3/status.json').read_text());status['v3_extraction_terminal_receipt']=extract
  limit_failure=extract.get('proxy_error')=='Transfer budget exceeded' or extract.get('bytes',0)>=extract.get('byte_limit',0)-16777216
  if not limit_failure:raise RuntimeError('V3 failed for another reason; no automatic restart: '+str(extract))
  assert not (N/'verified_HAL_v4_chain_20261009').exists(),'V4 already exists; duplicate forbidden'
  run('run_verified_HAL_v4_chain.py','recover_v4_with_existing_validated_cache')
  new=json.loads((N/'verified_HAL_v4_chain_20261009/status.json').read_text());status['v4_terminal_receipt']=new
  assert new['state']=='completed','V4 did not complete; no repeated recovery'
  run('review_remote_conservation_taxa_v4.py','review_v4_taxa');status.update(state='completed',version='v4',stage='validated_taxa_review')
 else:raise RuntimeError('Existing chain stopped without terminal receipt; investigate before restart')
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:status['utc_end']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
