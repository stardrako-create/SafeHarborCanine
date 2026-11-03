from pathlib import Path
import pysam,json,collections,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_SPADES_REVIEW_20261004';out={}
for label in ['ROS','UU']:
 hits=collections.defaultdict(list)
 with pysam.AlignmentFile(str(O/(label+'.sam')),'r') as f:
  for x in f:
   if x.is_unmapped:continue
   hits[x.query_name].append({'ref':x.reference_name,'start0':x.reference_start,'end0':x.reference_end,'cigar':x.cigarstring,'MAPQ':x.mapping_quality,'AS':x.get_tag('AS'),'NM':x.get_tag('NM'),'reverse':x.is_reverse,'secondary':x.is_secondary,'supplementary':x.is_supplementary})
 out[label]={'mapped_contig_names':len(hits),'junction_contig_hits':hits['NODE_20_length_175_cov_6.704082'],'primary_summary':{k:[x for x in v if not x['secondary'] and not x['supplementary']] for k,v in hits.items()}}
(O/'alignment_review.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'references':out,'limits':['BWA reported alignments are heuristic','Contig from same reads is not independent evidence','Partial repeated alignments do not establish distant origin or complete inserted allele']},indent=2));print(json.dumps({k:{'mapped':v['mapped_contig_names'],'junction_hits':v['junction_contig_hits']} for k,v in out.items()}))