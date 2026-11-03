from pathlib import Path
import http.server,threading,urllib.request,subprocess,json,time,traceback,re,os
N=Path(__file__).resolve().parent;O=N/'remote_HAL_verified_20261009_v2';integrity=json.loads((N/'cache_integrity_20261009/status.json').read_text());assert integrity['state']=='completed' and integrity['checked']==integrity['blocks'] and not integrity['errors'];O.mkdir(exist_ok=False)
URL='https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.hal';counter={'bytes':0,'requests':0,'cache_hits':0};disk=Path('/mnt/d/Jin2024_work/hal_range_cache_verified_20261009');disk.mkdir(exist_ok=True);cache={};lock=threading.Lock();status={'state':'running','source':URL,'byte_limit':268435456}
def durable(path,data):
 temporary=path.with_name(path.name+f'.tmp.{os.getpid()}.{threading.get_ident()}')
 with temporary.open('wb') as f:f.write(data);f.flush();os.fsync(f.fileno())
 for attempt in range(100):
  try:temporary.replace(path);return
  except PermissionError:
   if attempt==99:raise
   time.sleep(.05)
last_saved=0
def save(force=False):
 global last_saved
 now=time.monotonic()
 if not force and now-last_saved<2:return
 durable(O/'status.json',json.dumps(dict(status,**counter),indent=2).encode());last_saved=now
class Handler(http.server.BaseHTTPRequestHandler):
 protocol_version='HTTP/1.1'
 def log_message(self,*args):pass
 def do_HEAD(self):self.respond(True)
 def do_GET(self):self.respond(False)
 def respond(self,head):
  try:
   if self.path!='/alignment.hal':self.send_error(404);return
   headers={}
   if not head:
    r=self.headers.get('Range','');m=re.fullmatch(r'bytes=(\d+)-(\d+)',r)
    if not m: self.send_error(400,'Range required');return
    lo,hi=map(int,m.groups());size=hi-lo+1
    if size<1 or size>16777216: self.send_error(413);return
    with lock:
     if counter['bytes']+size>268435456:self.send_error(507);return
    headers['Range']=r
   if head:
    with urllib.request.urlopen(urllib.request.Request(URL,method='HEAD'),timeout=30) as resp:
     meta=dict(resp.headers)
    signature={k:meta.get(k) for k in ['Content-Length','Last-Modified','ETag']};signature['url']=URL
    manifest=disk/'source.json'
    if manifest.exists():assert json.loads(manifest.read_text())==signature,'Remote source changed'
    else:manifest.write_text(json.dumps(signature,indent=2))
    self.send_response(200)
    for k in ['Content-Length','Last-Modified','ETag']:
     if meta.get(k):self.send_header(k,meta[k])
    self.send_header('Accept-Ranges','bytes');self.end_headers();return
   total=864878803587;block=65536;parts=[]
   for begin in range((lo//block)*block,hi+1,block):
    end=min(total-1,begin+block-1)
    blockfile=disk/(str(begin)+'.bin')
    if begin not in cache and blockfile.exists():
     data=blockfile.read_bytes();assert len(data)==end-begin+1;cache[begin]=data;counter['cache_hits']+=1;save()
    if begin not in cache:
     with lock:
      if counter['bytes']+end-begin+1>268435456:raise RuntimeError('Transfer budget exceeded')
     req=urllib.request.Request(URL,headers={'Range':f'bytes={begin}-{end}'})
     with urllib.request.urlopen(req,timeout=30) as resp:
      if resp.status!=206 or resp.headers.get('Content-Range')!=f'bytes {begin}-{end}/{total}':raise RuntimeError('Invalid upstream range')
      assert resp.headers.get('ETag')==integrity['source']['ETag'];data=resp.read(end-begin+2)
     if len(data)!=end-begin+1:raise RuntimeError('Range size mismatch')
     cache[begin]=data;durable(blockfile,data)
     with lock:counter['bytes']+=len(data);counter['requests']+=1;save()
    parts.append(cache[begin][max(lo,begin)-begin:min(hi+1,end+1)-begin])
   data=b''.join(parts);assert len(data)==size
   self.send_response(206);self.send_header('Content-Length',str(size));self.send_header('Content-Range',f'bytes {lo}-{hi}/{total}');self.send_header('Accept-Ranges','bytes');self.end_headers();self.wfile.write(data)

  except Exception as e:
   status['proxy_error']=str(e);save();self.send_error(502,str(e))
server=http.server.HTTPServer(('127.0.0.1',0),Handler);threading.Thread(target=server.serve_forever,daemon=True).start()
mount=Path('/home/stardrako/safeharbor_hal_verified_20261009_v2');mount.mkdir(exist_ok=True);proc=None
try:
 log=(O/'httpfs.log').open('w');cmd=['/mnt/d/Jin2024_work/httpfs_probe_20261009/extracted/usr/bin/httpfs2','-f','-t','40',f'http://127.0.0.1:{server.server_port}/alignment.hal',str(mount),'-o','ro'];(O/'command.json').write_text(json.dumps(cmd));proc=subprocess.Popen(cmd,stdout=log,stderr=log)
 for _ in range(20):
  if os.path.ismount(mount):break
  if proc.poll() is not None:raise RuntimeError('httpfs exited '+str(proc.returncode))
  time.sleep(.5)
 assert os.path.ismount(mount),'Mount did not become ready'
 targets=O/'targets.bed';targets.write_text('chr27\t18845873\t18846873\nchr3\t2199848\t2200848\nchr10\t42483814\t42484814\n')
 cmd=['/home/stardrako/miniforge3/envs/cactus/bin/hal2maf','--refGenome','Canis_lupus_familiaris','--refTargets',str(targets),'--noAncestors','--noDupes',str(mount/'alignment.hal'),str(O/'three_regions.maf')];(O/'extraction_command.json').write_text(json.dumps(cmd))
 with (O/'halStats.out').open('w') as out,(O/'halStats.err').open('w') as err:r=subprocess.run(cmd,stdout=out,stderr=err,timeout=7200)
 assert r.returncode==0,'hal2maf failed';status['state']='completed'
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:
 if os.path.ismount(mount):subprocess.run(['fusermount','-u',str(mount)],capture_output=True,timeout=10)
 if proc and proc.poll() is None:proc.terminate()
 server.shutdown();save(True);print(json.dumps(dict(status,**counter)))


