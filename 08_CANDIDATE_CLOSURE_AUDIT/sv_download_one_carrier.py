"""Download only PBGV000010 mates, sequentially; preserve all partial/unvalidated data."""
from pathlib import Path
import csv,json,hashlib,urllib.request,shutil,time,os,msvcrt
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT';D=O/'SV_PBGV000010';D.mkdir(exist_ok=True)
lock=(D/'worker.lock').open('a+b');lock.seek(0)
try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
except OSError:raise SystemExit('Another worker holds the download lock')
state={'sample':'PBGV000010','run':'SRR15734832','pid':os.getpid(),'status':'starting','validated_files':[]}
def save(**values):
 state.update(values);state['updated_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
 tmp=D/'status.tmp';tmp.write_text(json.dumps(state,indent=2));os.replace(tmp,D/'status.json')
def md5(path):
 h=hashlib.md5()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
try:
 rows=list(csv.DictReader((O/'SRR15734832_ENA_verified_manifest.tsv').open(),delimiter='\t'));assert len(rows)==1
 r=rows[0];assert r['sample_accession']=='SAMN21036425' and r['run_accession']=='SRR15734832' and r['library_strategy']=='WGS' and r['library_layout']=='PAIRED'
 entries=list(zip(r['fastq_ftp'].split(';'),map(int,r['fastq_bytes'].split(';')),r['fastq_md5'].split(';')));assert len(entries)==2
 for remote,size,expected in entries:
  target=D/remote.rsplit('/',1)[1];part=target.with_suffix(target.suffix+'.part')
  if target.exists():
   assert target.stat().st_size==size and md5(target)==expected,'Existing final file does not validate; preserved'
  else:
   offset=part.stat().st_size if part.exists() else 0;assert offset<=size
   assert shutil.disk_usage(D).free>(size-offset)+45*1024**3,'Insufficient disk reserve for this mate'
   if offset<size:
    req=urllib.request.Request('https://'+remote,headers={'Range':f'bytes={offset}-'} if offset else {})
    with urllib.request.urlopen(req,timeout=90) as response:
     if offset:assert response.status==206 and response.headers['Content-Range'].startswith(f'bytes {offset}-'),'Server did not honor resume; partial preserved'
     with part.open('ab' if offset else 'wb') as f:
      last=0
      while True:
       b=response.read(4*1024*1024)
       if not b:break
       if shutil.disk_usage(D).free<45*1024**3:raise RuntimeError('Disk reserve reached; partial preserved')
       f.write(b);offset+=len(b)
       if time.monotonic()-last>15:
        save(status='downloading',file=target.name,bytes_received=offset,expected_bytes=size,free_bytes=shutil.disk_usage(D).free);last=time.monotonic()
   save(status='validating_md5',file=target.name,bytes_received=part.stat().st_size,expected_bytes=size)
   assert part.stat().st_size==size,'Wrong length; partial preserved'
   actual=md5(part);assert actual==expected,'MD5 mismatch; partial preserved'
   part.rename(target)
  receipt={'file':target.name,'bytes':size,'md5':expected,'validated':True}
  (D/(target.name+'.validated.json')).write_text(json.dumps(receipt,indent=2));state['validated_files'].append(receipt)
 save(status='download_complete',note='Both ENA sizes and MD5 verified. No alignment or biological validation performed.')
except Exception as e:
 save(status='failed',error=repr(e));raise
finally:
 lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1);lock.close()
