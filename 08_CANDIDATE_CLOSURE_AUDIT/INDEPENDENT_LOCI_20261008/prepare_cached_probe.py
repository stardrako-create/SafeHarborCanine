from pathlib import Path
N=Path(__file__).resolve().parent
s=(N/'probe_remote_HAL_v2.py').read_text(encoding='utf-8-sig').replace("remote_HAL_probe_20261009_v2'","remote_HAL_cached_20261009'").replace('safeharbor_hal_remote_20261009','safeharbor_hal_cached_20261009')
s=s.replace("lock=threading.Lock();status=","cache={};lock=threading.Lock();status=")
a=s.index('   req=urllib.request.Request');b=s.index('\n  except Exception as e:',a)
s=s[:a]+'''   if head:
    with urllib.request.urlopen(urllib.request.Request(URL,method='HEAD'),timeout=30) as resp:
     meta=dict(resp.headers)
    self.send_response(200)
    for k in ['Content-Length','Last-Modified','ETag']:
     if meta.get(k):self.send_header(k,meta[k])
    self.send_header('Accept-Ranges','bytes');self.end_headers();return
   total=864878803587;block=65536;parts=[]
   for begin in range((lo//block)*block,hi+1,block):
    end=min(total-1,begin+block-1)
    if begin not in cache:
     with lock:
      if counter['bytes']+end-begin+1>268435456:raise RuntimeError('Transfer budget exceeded')
     req=urllib.request.Request(URL,headers={'Range':f'bytes={begin}-{end}'})
     with urllib.request.urlopen(req,timeout=30) as resp:
      if resp.status!=206 or resp.headers.get('Content-Range')!=f'bytes {begin}-{end}/{total}':raise RuntimeError('Invalid upstream range')
      data=resp.read(end-begin+2)
     if len(data)!=end-begin+1:raise RuntimeError('Range size mismatch')
     cache[begin]=data
     with lock:counter['bytes']+=len(data);counter['requests']+=1;save()
    parts.append(cache[begin][max(lo,begin)-begin:min(hi+1,end+1)-begin])
   data=b''.join(parts);assert len(data)==size
   self.send_response(206);self.send_header('Content-Length',str(size));self.send_header('Content-Range',f'bytes {lo}-{hi}/{total}');self.send_header('Accept-Ranges','bytes');self.end_headers();self.wfile.write(data)
''' + s[b:]
s=s.replace("cmd=['/home/stardrako/miniforge3/envs/cactus/bin/halStats','--genomes',str(mount/'alignment.hal')]", "targets=O/'targets.bed';targets.write_text('chr27\\t18845873\\t18846873\\nchr3\\t2199848\\t2200848\\nchr10\\t42483814\\t42484814\\n')\n cmd=['/home/stardrako/miniforge3/envs/cactus/bin/hal2maf','--refGenome','Canis_lupus_familiaris','--refTargets',str(targets),'--noAncestors','--noDupes',str(mount/'alignment.hal'),str(O/'three_regions.maf')];(O/'extraction_command.json').write_text(json.dumps(cmd))")
s=s.replace('timeout=120','timeout=3600').replace("'halStats failed'","'hal2maf failed'")
(N/'extract_remote_HAL_cached.py').write_text(s)
