from pathlib import Path
import pysam,json,collections,hashlib,datetime
A=Path(__file__).resolve().parent;I=A/'SV_PBGV000012/ROS_expanded_junction_20261001';O=A/'W11_BRANCH_PHASE_20261002';O.mkdir(exist_ok=False)
r=json.loads((I/'review.json').read_text());keys={x['read'] for x in r['rows']};rows=[];counts={p:collections.Counter() for p in [1062,1073]}
with pysam.AlignmentFile(str(I/'joint.sam'),'r') as f:
 for x in f:
  if x.query_name not in keys or x.reference_name!='junction_L4707939_R33613308':continue
  assert x.cigarstring=='150M'
  bases={}
  for q,p in x.get_aligned_pairs(matches_only=True):
   if p in counts:
    b=x.query_sequence[q];qual=x.query_qualities[q];bases[str(p)]={'base':b,'BQ':qual}
    if qual>=20:counts[p][b]+=1
  rows.append({'read':x.query_name,'start0':x.reference_start,'end0':x.reference_end,'branch_bases':bases})
assert len(rows)==12
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'coordinate_system':'0-based synthetic junction, right starts1000','BQ20_counts':{str(k):dict(v) for k,v in counts.items()},'rows':rows,'limits':['Only12 previously selected spanning reads; same evidence reused','Base at branch is linked within one read to junction, not to distant genomic locus','Does not extend physical linkage beyond maximum122bp right anchor','No inference of independent molecules or structural genotype']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),I/'review.json',I/'joint.sam']},indent=2));print(json.dumps(out['BQ20_counts']))