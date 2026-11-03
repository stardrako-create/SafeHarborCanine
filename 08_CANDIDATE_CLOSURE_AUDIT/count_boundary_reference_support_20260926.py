from pathlib import Path
import pysam,json,datetime,collections
A=Path(__file__).resolve().parent;P=A/'SV_PBGV000010/alignment_chunked';O=A/'SV_CALL_AUDIT_20260926';rows=[]
with pysam.AlignmentFile(str(P/'selected_markdup.bam'),'rb') as f:
 for label,pos in [('left',4707939),('right',33613308)]:
  span=set();clip=set();eligible=set()
  for r in f.fetch('NC_049253.1',pos-200,pos+200):
   if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality<30:continue
   eligible.add(r.query_name)
   if any(b<=pos-20 and e>=pos+20 for b,e in r.get_blocks()):span.add(r.query_name)
   if r.cigartuples and ((r.reference_end==pos and r.cigartuples[-1][0]==4 and r.cigartuples[-1][1]>=20) or (r.reference_start==pos and r.cigartuples[0][0]==4 and r.cigartuples[0][1]>=20)):clip.add(r.query_name)
  rows.append(dict(boundary=label,pos0=pos,eligible_names_within_200bp=len(eligible),reference_spanning_names_20bp_each_side=len(span),exact_boundary_softclip_names_ge20bp=len(clip),both_categories=len(span&clip),reference_names=sorted(span),softclip_names=sorted(clip)))
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rows':rows,'limits':['Read-name groups not verified independent molecules','Reference span criterion tests CIGAR block continuity, not perfect base identity','Selected chromosome BAM; counts are descriptive, not genotype likelihood or variant allele fraction','Soft clips are not automatically support for the distant junction']};(O/'boundary_reference_support.json').write_text(json.dumps(r,indent=2));print(json.dumps([{k:v for k,v in x.items() if not k.endswith('_names')} for x in rows]))
