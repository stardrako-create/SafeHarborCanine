from pathlib import Path
import urllib.request,hashlib,json
O=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino\06_v1.21.6')
u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244116/suppl/GSE244116_filtered_feature_bc_matrix.h5'
p=O/'GSE244116_filtered_feature_bc_matrix.h5';tmp=p.with_suffix('.partial');n=0;sha=hashlib.sha256()
with urllib.request.urlopen(u,timeout=45) as inp,tmp.open('wb') as out:
 expected=int(inp.headers['Content-Length'])
 while True:
  b=inp.read(1024*1024)
  if not b:break
  out.write(b);sha.update(b);n+=len(b)
assert n==expected;tmp.replace(p)
(O/'GSE244116_download.json').write_text(json.dumps(dict(url=u,bytes=n,sha256=sha.hexdigest(),assembly='GCF_000002285.5 canFam6',donors=1,context='primary treatment-naive osteosarcoma; immune nuclei are tumor-infiltrating',fragments_url='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7807nnn/GSM7807442/suppl/GSM7807442_atac_fragments.tsv.gz',fragments_downloaded=False),indent=2))
print('Matrix verified',n,sha.hexdigest())
