from pathlib import Path
import json,pysam,hashlib
O=Path(__file__).parent
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
rows=json.loads((O/'UU_projection/projection_review.json').read_text())['results']
results=[]
with pysam.FastaFile(str(O/'UU_projection/controls_ROS.fa')) as ros,pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa') as uu:
 for w in rows:
  h=next(x for x in w['hits'] if x['primary'])
  if h['cigar']!='1000M':continue
  a=ros.fetch(w['id']).upper();b=uu.fetch(h['chrom'],h['start0'],h['end0']).upper()
  aligned=rc(b) if h['strand']=='-' else b
  diffs=[]
  for i,(x,y) in enumerate(zip(a,aligned)):
   if x!=y:diffs.append(dict(ROS_offset0=i,ROS_base=x,UU_base_in_ROS_orientation=y,UU_position0=h['end0']-1-i if h['strand']=='-' else h['start0']+i))
  assert len(a)==len(aligned)==1000 and len(diffs)==h['NM']
  results.append(dict(id=w['id'],mapping=h,differences=diffs,ROS_sequence_sha256=hashlib.sha256(a.encode()).hexdigest(),UU_oriented_sequence_sha256=hashlib.sha256(aligned.encode()).hexdigest()))
(O/'ungapped_reference_differences.json').write_text(json.dumps(dict(results=results,limitations=['Assembly differences are not individual genotypes','Nonidentical REF must be reconciled per allele before projection','Gapped alignments require separate CIGAR mapping']),indent=2))
print(json.dumps({'windows':len(results),'mismatched_windows':sum(bool(x['differences']) for x in results)}))