from pathlib import Path
import pysam,json,datetime,collections
A=Path(__file__).resolve().parent;O=A/'W11_RAW_RECRUIT_20261005';names={'SRR14750332.53211536','SRR14750332.171739498'};hits=collections.defaultdict(list)
with pysam.AlignmentFile(str(O/'joint_pairs.sam'),'r') as f:
 for r in f:
  if r.query_name not in names or r.is_unmapped:continue
  key=r.query_name+('_mate1' if r.is_read1 else '_mate2');h={'ref':r.reference_name,'start0':r.reference_start,'end0':r.reference_end,'cigar':r.cigarstring,'MAPQ':r.mapping_quality,'AS':r.get_tag('AS'),'NM':r.get_tag('NM'),'reverse':r.is_reverse,'secondary':r.is_secondary,'supplementary':r.is_supplementary}
  if r.reference_name.startswith('junction_') and any(s<=980 and e>=1020 for s,e in r.get_blocks()):
   qs=[q for q,p in r.get_aligned_pairs(matches_only=True) if 980<=p<1020];h['spans_20_each']=True;h['anchor40_minBQ']=min(r.query_qualities[q] for q in qs) if r.query_qualities else None
  hits[key].append(h)
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'hits':dict(hits),'limits':['Two newly retrieved pairs, not yet molecule-deduplicated with full dataset','Small selected paired alignment not calibrated insert-size model','MAPQ affected by near-identical synthetic junction alternatives','Junction support does not establish distal source']};(O/'new_pairs_review.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:[h for h in v if h.get('spans_20_each') or not(h['secondary'] or h['supplementary'])] for k,v in hits.items()}))