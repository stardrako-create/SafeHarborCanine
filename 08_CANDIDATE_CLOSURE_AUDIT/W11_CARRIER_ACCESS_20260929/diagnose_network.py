import pathlib,subprocess,sys,concurrent.futures,json,datetime
hosts=['sra-pub-run-odp.s3.amazonaws.com','locate.ncbi.nlm.nih.gov','www.ebi.ac.uk','kiddlabshare.med.umich.edu','sra-downloadb.be-md.ncbi.nlm.nih.gov']
code="""import socket,json,sys,time
h=sys.argv[1];t=time.time();a=socket.getaddrinfo(h,443,type=socket.SOCK_STREAM);o=[]
for fam,typ,pro,_,addr in a[:3]:
 s=socket.socket(fam,typ,pro);s.settimeout(2)
 try:s.connect(addr);v='connected'
 except Exception as e:v=str(e)
 finally:s.close()
 o.append({'address':addr[0],'family':int(fam),'tcp':v})
print(json.dumps({'host':h,'elapsed':time.time()-t,'results':o}))
"""
def f(h):
 try:
  r=subprocess.run([sys.executable,'-c',code,h],capture_output=True,text=True,timeout=12);return {'host':h,'out':r.stdout,'error':r.stderr}
 except subprocess.TimeoutExpired:return {'host':h,'error':'DNS/connect hard timeout 12s'}
with concurrent.futures.ThreadPoolExecutor(5) as x:r=list(x.map(f,hosts))
p=pathlib.Path(__file__).parent/('network_diagnosis_'+datetime.datetime.now(datetime.timezone.utc).strftime('%H%M%S')+'.json');p.write_text(json.dumps(r,indent=2));print(json.dumps(r))
