from pathlib import Path
import json,urllib.request,concurrent.futures,hashlib,os,datetime,traceback
N=Path(__file__).resolve().parent;O=N/'cache_integrity_20261009';O.mkdir(exist_ok=False)
old=Path('/mnt/d/Jin2024_work/hal_range_cache_20261009');new=Path('/mnt/d/Jin2024_work/hal_range_cache_verified_20261009');new.mkdir(exist_ok=False)
source=json.loads((old/'source.json').read_text());url=source['url'];total=int(source['Content-Length'])
with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=30) as r:
 assert all(r.headers.get(k)==source[k] for k in ['Content-Length','ETag','Last-Modified'])
(new/'source.json').write_text(json.dumps(source,indent=2))
files=sorted(old.glob('*.bin'),key=lambda p:int(p.stem));groups=[]
for p in files:
 s=int(p.stem);assert s%65536==0
 if groups and s==groups[-1][-1]+65536 and len(groups[-1])<16:groups[-1].append(s)
 else:groups.append([s])
status=dict(state='running',source=source,blocks=len(files),groups=len(groups),checked=0,mismatched=0,errors=[],bytes=0)
def save():
 tmp=O/'status.tmp'
 with tmp.open('w') as f:json.dump(status,f,indent=2);f.flush();os.fsync(f.fileno())
 tmp.replace(O/'status.json')
save()
def check(group):
 lo=group[0];hi=min(total-1,group[-1]+65535);req=urllib.request.Request(url,headers={'Range':f'bytes={lo}-{hi}'})
 with urllib.request.urlopen(req,timeout=60) as r:
  assert r.status==206 and r.headers.get('Content-Range')==f'bytes {lo}-{hi}/{total}'
  assert r.headers.get('ETag')==source['ETag'];data=r.read(hi-lo+2)
 assert len(data)==hi-lo+1
 results=[]
 for s in group:
  block=data[s-lo:min(s-lo+65536,len(data))];before=(old/(str(s)+'.bin')).read_bytes();path=new/(str(s)+'.bin');tmp=path.with_suffix('.tmp')
  with tmp.open('wb') as f:f.write(block);f.flush();os.fsync(f.fileno())
  tmp.replace(path);results.append(dict(offset=s,size=len(block),matched=before==block,old_sha256=hashlib.sha256(before).hexdigest(),verified_sha256=hashlib.sha256(block).hexdigest()))
 return results,len(data)
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 tasks={pool.submit(check,g):g for g in groups}
 for task in concurrent.futures.as_completed(tasks):
  try:
   rr,size=task.result();results.extend(rr);status['checked']+=len(rr);status['mismatched']+=sum(not r['matched'] for r in rr);status['bytes']+=size
  except Exception as e:status['errors'].append(dict(offset=tasks[task][0],error=str(e)))
  save()
(O/'block_checks.json').write_text(json.dumps(sorted(results,key=lambda r:r['offset']),indent=2));status['state']='completed' if not status['errors'] and status['checked']==status['blocks'] else 'failed';save();print(json.dumps(status))
