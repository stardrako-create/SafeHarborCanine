from pathlib import Path
import requests,json
N=Path(__file__).resolve().parent;O=N/'local_gene_ortholog_probe_20261009';O.mkdir(exist_ok=False);url='https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/486731/orthologs';r=requests.get(url,params={'taxon_filters':'9606','page_size':1000},timeout=60);(O/'raw_response.txt').write_bytes(r.content);(O/'request.json').write_text(json.dumps(dict(url=r.url,status_code=r.status_code,headers=dict(r.headers)),indent=2));print(r.status_code);print(r.text[:2500])
