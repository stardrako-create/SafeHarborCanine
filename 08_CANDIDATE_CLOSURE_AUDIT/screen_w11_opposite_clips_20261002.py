from pathlib import Path
import pysam,json,collections,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_OPPOSITE_CLIPS_20261002';O.mkdir(exist_ok=False)
B=A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam';rows=[]
with pysam.AlignmentFile(str(B),'rb') as f:
 for r in f.fetch('NC_051836.1',35341619,35342633):
  if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality<30:continue
  for side,op in [('left',r.cigartuples[0]),('right',r.cigartuples[-1])]:
   if op[0]!=4 or op[1]<20:continue
   n=op[1];seq=r.query_sequence[:n] if side=='left' else r.query_sequence[-n:];q=list(r.query_qualities[:n] if side=='left' else r.query_qualities[-n:]);boundary=r.reference_start if side=='left' else r.reference_end
   rows.append({'name':r.query_name,'mate':1 if r.is_read1 else 2,'side':side,'boundary0':boundary,'clip':seq,'clip_BQ20_bases':sum(x>=20 for x in q),'clip_all_BQ20':min(q)>=20,'MAPQ':r.mapping_quality,'cigar':r.cigarstring,'reverse':r.is_reverse,'start0':r.reference_start,'end0':r.reference_end})
groups=[]
for side,b in sorted({(x['side'],x['boundary0']) for x in rows}):
 rr=[x for x in rows if x['side']==side and x['boundary0']==b];groups.append({'side':side,'boundary0':b,'reads':len(rr),'distinct_names':len({x['name'] for x in rr}),'all_clip_BQ20':sum(x['clip_all_BQ20'] for x in rr),'max_clip':max(len(x['clip']) for x in rr)})
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'groups':groups,'rows':rows,'limits':['Local primary MAPQ30 nonduplicate QCpass reads; softclip>=20','Clips alone do not establish insertion or distant origin','Coordinates are alignment boundaries, not validated breakpoints']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'script.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest());print(json.dumps(groups))