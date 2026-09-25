from pathlib import Path
import urllib.request,hashlib,json,time
O=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino\06_v1.21.7');O.mkdir(exist_ok=True)
urls=['https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToCanFam6.over.chain.gz','https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/liftOver/canFam6ToGCF_014441545.1.over.chain.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7807nnn/GSM7807442/suppl/GSM7807442_atac_fragments.tsv.gz']
for u in urls:
 p=O/u.rsplit('/',1)[1];tmp=p.with_suffix(p.suffix+'.partial');t=time.time();n=0;sha=hashlib.sha256()
 if p.exists():continue
 with urllib.request.urlopen(u,timeout=60) as inp,tmp.open('wb') as out:
  expected=int(inp.headers['Content-Length']);last=0
  while True:
   b=inp.read(1024*1024)
   if not b:break
   out.write(b);sha.update(b);n+=len(b)
   if time.time()-last>30:print(p.name,n,expected,flush=True);last=time.time()
 assert n==expected;tmp.replace(p)
 (O/(p.name+'.download.json')).write_text(json.dumps(dict(url=u,bytes=n,sha256=sha.hexdigest(),seconds=time.time()-t),indent=2))
 print('Complete',p.name,flush=True)
