from pathlib import Path
import http.server,threading,urllib.request,subprocess,json,time,traceback,re,os
N=Path(__file__).resolve().parent;O=N/'remote_HAL_probe_20261009_v2';O.mkdir(exist_ok=False)
URL='https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2.hal';counter={'bytes':0,'requests':0};lock=threading.Lock();status={'state':'running','source':URL,'byte_limit':268435456}
def save(): (O/'status.json').write_text(json.dumps(dict(status,**counter),indent=2))
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
   req=urllib.request.Request(URL,headers=headers,method='HEAD' if head else 'GET')
   with urllib.request.urlopen(req,timeout=30) as resp:
    if not head and resp.status!=206:raise RuntimeError('Upstream ignored Range')
    data=b'' if head else resp.read(size+1)
    if not head and len(data)!=size:raise RuntimeError('Range size mismatch')
    self.send_response(resp.status)
    for k in ['Content-Length','Content-Range','Last-Modified','ETag']:
     if resp.headers.get(k):self.send_header(k,resp.headers[k])
    self.send_header('Accept-Ranges','bytes');self.end_headers()
    if data:self.wfile.write(data)
    with lock:counter['bytes']+=len(data);counter['requests']+=1;save()
  except Exception as e:
   status['proxy_error']=str(e);save();self.send_error(502,str(e))
server=http.server.ThreadingHTTPServer(('127.0.0.1',0),Handler);threading.Thread(target=server.serve_forever,daemon=True).start()
mount=Path('/home/stardrako/safeharbor_hal_remote_20261009');mount.mkdir(exist_ok=True);proc=None
try:
 log=(O/'httpfs.log').open('w');cmd=['/mnt/d/Jin2024_work/httpfs_probe_20261009/extracted/usr/bin/httpfs2','-f','-t','40',f'http://127.0.0.1:{server.server_port}/alignment.hal',str(mount),'-o','ro'];(O/'command.json').write_text(json.dumps(cmd));proc=subprocess.Popen(cmd,stdout=log,stderr=log)
 for _ in range(20):
  if os.path.ismount(mount):break
  if proc.poll() is not None:raise RuntimeError('httpfs exited '+str(proc.returncode))
  time.sleep(.5)
 assert os.path.ismount(mount),'Mount did not become ready'
 cmd=['/home/stardrako/miniforge3/envs/cactus/bin/halStats','--genomes',str(mount/'alignment.hal')]
 with (O/'halStats.out').open('w') as out,(O/'halStats.err').open('w') as err:r=subprocess.run(cmd,stdout=out,stderr=err,timeout=120)
 assert r.returncode==0,'halStats failed';status['state']='completed'
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:
 if os.path.ismount(mount):subprocess.run(['fusermount','-u',str(mount)],capture_output=True,timeout=10)
 if proc and proc.poll() is None:proc.terminate()
 server.shutdown();save();print(json.dumps(dict(status,**counter)))

