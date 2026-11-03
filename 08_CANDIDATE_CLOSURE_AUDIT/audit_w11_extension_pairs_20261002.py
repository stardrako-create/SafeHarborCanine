from pathlib import Path
import json,pysam,datetime,hashlib
A=Path(__file__).resolve().parent;I=A/'W11_OVERLAP_20261002/review.json';O=A/'W11_EXTENSION_PAIRS_20261002';O.mkdir(exist_ok=False)
r=json.loads(I.read_text()); candidates={c['read'] for row in r['rows'] for c in row['candidates'] if c['right_end_if_correct']>122}; names={x.rsplit('_mate',1)[0] for x in candidates};records=[]
B=A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam'
with pysam.AlignmentFile(str(B),'rb') as f:
 for x in f:
  if x.query_name not in names or x.is_secondary or x.is_supplementary:continue
  records.append(dict(name=x.query_name,mate=1 if x.is_read1 else 2,ref=x.reference_name,start0=x.reference_start,end0=x.reference_end,cigar=x.cigarstring,mapq=x.mapping_quality,reverse=x.is_reverse,duplicate=x.is_duplicate,qcfail=x.is_qcfail,flag=x.flag,SA=x.get_tag('SA') if x.has_tag('SA') else None))
assert len(candidates)==5
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'candidate_reads':sorted(candidates),'primary_records':records,'limits':['Retained chr32-linked BAM only','MAPQ and orientation are observations, not calibrated pair likelihood','A mate near the breakpoint can anchor a repetitive read without assigning its distal source','Single-step overlap candidates are not assembled molecules']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),I]},indent=2));print(json.dumps(out))