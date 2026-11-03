import pathlib,urllib.request,json,datetime,shutil,hashlib,msvcrt,os
ROOT=pathlib.Path(__file__).parent
DEST=pathlib.Path.home()/'SafeHarborData'/'SRR14750332';DEST.mkdir(parents=True,exist_ok=True)
SIZE=26235010732; MD5='035a1a8a3add18107a9cc8763acb6048';RESERVE=45*1024**3
URL='https://sra-pub-run-odp.s3.dualstack.us-east-1.amazonaws.com/sra/SRR14750332/SRR14750332'
part=DEST/'SRR14750332.sra.part';final=DEST/'SRR14750332.sra'
lock=open(DEST/'download.lock','a+b');lock.seek(0);lock.write(b'0');lock.flush();lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
def status(state,**kw):
 d={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'state':state,'url':URL,'destination':str(final),'expected_bytes':SIZE,'expected_md5':MD5,**kw};tmp=ROOT/'download_status.tmp';tmp.write_text(json.dumps(d,indent=2));tmp.replace(ROOT/'download_status.json')
try:
 if final.exists():raise RuntimeError('Final already exists; preserve and verify separately')
 offset=part.stat().st_size if part.exists() else 0
 if shutil.disk_usage(DEST).free < SIZE-offset+RESERVE:raise RuntimeError('Insufficient space for remaining data plus 45 GiB reserve')
 status('downloading',bytes=offset)
 with part.open('ab') as f:
  while offset<SIZE:
   if shutil.disk_usage(DEST).free<RESERVE+8*1024**2:raise RuntimeError('Reserve reached')
   end=min(offset+8*1024**2,SIZE)-1
   req=urllib.request.Request(URL,headers={'Range':f'bytes={offset}-{end}'})
   with urllib.request.urlopen(req,timeout=30) as r:
    if r.status!=206 or r.headers.get('Content-Range')!=f'bytes {offset}-{end}/{SIZE}':raise RuntimeError('Unexpected range response')
    b=r.read(end-offset+2)
   if len(b)!=end-offset+1:raise RuntimeError('Short range; no chunk appended')
   f.write(b);f.flush();offset+=len(b);status('downloading',bytes=offset)
 status('verifying_md5',bytes=offset)
 h=hashlib.md5()
 with part.open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 if part.stat().st_size!=SIZE or h.hexdigest()!=MD5:raise RuntimeError('Size/MD5 mismatch; preserve part')
 part.rename(final);status('download_verified',bytes=SIZE,observed_md5=h.hexdigest())
except Exception as e:
 status('stopped',error=str(e),bytes=part.stat().st_size if part.exists() else 0);raise
finally:lock.close()
