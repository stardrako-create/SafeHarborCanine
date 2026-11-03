import concurrent.futures,subprocess,sys,json,datetime,pathlib,shutil
ROOT=pathlib.Path(__file__).parent
URLS={
 'original_sra_s3':'https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR12330331/SRR12330331',
 'ncbi_resolver':'https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=SRR12330331&location=ncbi',
 'ena':'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR12330331&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5&format=tsv'}
CHILD='''import urllib.request,json,sys
u=sys.argv[1]
req=urllib.request.Request(u,headers={'Range':'bytes=0-1023'})
with urllib.request.urlopen(req,timeout=8) as r:
 b=r.read(65536 if 'retrieve?' in u or 'filereport?' in u else 1024)
 print(json.dumps({'status':r.status,'content_range':r.headers.get('Content-Range'),'bytes_read':len(b),'body':b.decode(errors='replace') if 'retrieve?' in u or 'filereport?' in u else None}))
'''
def probe(k,u):
 try:
  r=subprocess.run([sys.executable,'-c',CHILD,u],capture_output=True,text=True,timeout=18)
  return {'name':k,'url':u,'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr[-1500:]}
 except subprocess.TimeoutExpired:return {'name':k,'url':u,'status':'hard_timeout_18s_child_terminated'}
if __name__=='__main__':
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(lambda x:probe(*x),URLS.items()))
 report={'utc':stamp,'results':rows,'free_bytes':{d:shutil.disk_usage(d+':/').free for d in ['C','D']},'scientific_status':'no new reads validated; metadata access is not sequence validation'}
 out=ROOT/('bounded_access_GBGV_'+stamp+'.json');out.write_text(json.dumps(report,indent=2),encoding='utf-8')
 print(str(out));print(json.dumps(report))
