from pathlib import Path
import csv,json,urllib.request,concurrent.futures,subprocess,hashlib,time,threading,shutil
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.6';D=O/'raw';D.mkdir(parents=True,exist_ok=True)
manifest=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'));files=[];sources=[]
for r in manifest:
 url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession='+r['SRR']+'&result=read_run&fields=run_accession,experiment_accession,sample_accession,fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count,library_strategy,library_layout&format=tsv'
 with urllib.request.urlopen(url,timeout=40) as f:data=f.read()
 rows=list(csv.DictReader(data.decode().splitlines(),delimiter='\t'));assert len(rows)==1;x=rows[0];assert x['run_accession']==r['SRR'] and x['experiment_accession']==r['SRX'] and x['library_layout']=='PAIRED' and x['library_strategy']=='RNA-Seq'
 sources.append(dict(source=url,**x));urls=x['fastq_ftp'].split(';');sizes=x['fastq_bytes'].split(';');md5=x['fastq_md5'].split(';');assert len(urls)==len(sizes)==len(md5)==2
 for i,(u,n,h) in enumerate(zip(urls,sizes,md5),1):files.append(dict(sample=r['column'],run=r['SRR'],mate=i,url='https://'+u,bytes=int(n),md5=h,path=str(D/Path(u).name)))
(O/'ENA_verified_download_manifest.json').write_text(json.dumps(dict(files=files,metadata=sources),indent=2));assert shutil.disk_usage(D).free>sum(max(0,f['bytes']-(Path(f['path']).stat().st_size if Path(f['path']).exists() else 0)) for f in files)+50_000_000_000
lock=threading.Lock();state={f['run']+'_'+str(f['mate']):dict(status='pending',expected_bytes=f['bytes']) for f in files}
def save():
 temp=O/'download_status.tmp';temp.write_text(json.dumps(dict(updated=time.time(),files=state),indent=2));temp.replace(O/'download_status.json')
def get(f):
 key=f['run']+'_'+str(f['mate']);path=Path(f['path']);done=path.with_suffix(path.suffix+'.verified.json')
 if done.exists():
  previous=json.loads(done.read_text())
  if path.exists() and path.stat().st_size==f['bytes'] and previous.get('expected_md5')==f['md5']:
   with lock:state[key].update(status='verified',bytes=path.stat().st_size);save()
   return
 with lock:state[key]['status']='downloading';save()
 log=O/(key+'.download.log')
 with log.open('a') as out:
  for attempt in range(1,33):
   if path.exists() and path.stat().st_size==f['bytes']:break
   cmd=['curl.exe','-L','--fail','--retry','5','--retry-all-errors','--retry-delay','5','--connect-timeout','30','--max-time','1800','--continue-at','-','--output',str(path),f['url']]
   r=subprocess.run(cmd,stdout=out,stderr=out)
   with lock:state[key].update(attempt=attempt,curl_returncode=r.returncode,bytes=path.stat().st_size if path.exists() else 0);save()
   if r.returncode==0:break
  assert path.stat().st_size==f['bytes'],key+' size mismatch'
 with lock:state[key]['status']='checking_md5';save()
 h=hashlib.md5()
 with path.open('rb') as inp:
  for b in iter(lambda:inp.read(8*1024*1024),b''):h.update(b)
 assert h.hexdigest()==f['md5'],key+' MD5 mismatch'
 done.write_text(json.dumps(dict(expected_md5=f['md5'],observed_md5=h.hexdigest(),bytes=path.stat().st_size,verified_at=time.time()),indent=2))
 with lock:state[key].update(status='verified',bytes=path.stat().st_size);save()
 print(key,'verified',flush=True)
def guarded(f):
 try:get(f)
 except Exception as e:
  with lock:state[f['run']+'_'+str(f['mate'])].update(status='failed',error=str(e));save()
  raise
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(guarded,files))
print('All18 FASTQ files verified',flush=True)
