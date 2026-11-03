from pathlib import Path
import urllib.request,json,hashlib,datetime,os
N=Path(__file__).resolve().parent;O=N/'canFam3_reference_identity_20261009';O.mkdir(exist_ok=False)
results=[]
for label,chrom,s,e in [('1059','chr27',18845873,18846873),('1266','chr3',2199848,2200848),('1900','chr10',42483814,42484814)]:
 url=f'https://api.genome.ucsc.edu/getData/sequence?genome=canFam3;chrom={chrom};start={s};end={e}'
 try:
  with urllib.request.urlopen(url,timeout=90) as r:raw=r.read()
  (O/(label+'.response.json')).write_bytes(raw);obj=json.loads(raw);dna=obj['dna'];assert len(dna)==e-s and obj['start']==s and obj['end']==e and obj['chrom']==chrom
  (O/(label+'.fa')).write_text(f'>{chrom}:{s}-{e}\n{dna}\n')
  results.append(dict(id=label,url=url,start0=s,end0=e,chrom=chrom,length=len(dna),sequence_sha256=hashlib.sha256(dna.upper().encode()).hexdigest(),response_sha256=hashlib.sha256(raw).hexdigest(),state='retrieved'))
 except Exception as x:results.append(dict(id=label,url=url,state='failed',error=str(x)))
(O/'receipt.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),results=results,limitation='Independent canFam3 DNA retrieved; comparison to complete HAL MAF still required.'),indent=2));print(json.dumps(results))