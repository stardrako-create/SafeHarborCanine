from pathlib import Path
import json,hashlib,collections,datetime
N=Path(__file__).resolve().parent
checks=[]
for folder in ['variant_catalogue_query_drop_samples','variant_catalogue_query_remaining3']:
 for p in sorted((N/folder).glob('*__*.json')):
  j=json.loads(p.read_text());assert hashlib.sha256(p.with_suffix('.vcf').read_bytes()).hexdigest()==j['raw_vcf_sha256'];assert all(a['reference_verified'] for a in j['alleles']);checks.append(str(p.relative_to(N)))
blocks={r['id']:r for r in json.loads((N/'basewise_projection.json').read_text())['results']}
rows=json.loads((N/'normalized_complete_20261008/review.json').read_text())['results'];out=[]
for r in rows:
 b=blocks[r['window']];offset=r['start0']-int(r['window'].rsplit('_',2)[1]);end=offset+len(r['ref'])
 contained=any(x['ROS_offset_start0']<=offset<end<=x['ROS_offset_end0'] for x in b['exact_blocks'])
 out.append(dict(window=r['window'],id=r['id'],normalized_REF_within_exact_block=contained,start0=r['start0'],ref=r['ref'],alt=r['alt']))
res=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),checked_receipts=len(checks),raw_hash_and_reported_REF_checks_pass=True,normalized_count=len(rows),outside_exact_block=sum(not r['normalized_REF_within_exact_block'] for r in out),results=out,limitations=['within_exact_window in query receipts means alignment-envelope containment only for nonexact mappings; exact blocks are checked here separately.','Haplotype equivalence on ROS does not independently establish UU allele equivalence after crossing a block boundary.'])
(N/'normalization_boundary_review.json').write_text(json.dumps(res,indent=2));print(json.dumps({k:v for k,v in res.items() if k!='results'}))
