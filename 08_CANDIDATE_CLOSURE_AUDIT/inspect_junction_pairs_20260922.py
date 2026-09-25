from pathlib import Path
import pysam,json
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked';O=A/'junction_competition_20260922_v2';s=json.loads((O/'summary.json').read_text());names={r['read'].split('_mate')[0] for r in s['best_spanning_per_read']};rows=[]
with pysam.AlignmentFile(str(A/'selected_markdup.bam'),'rb') as f:
 for r in f.fetch(until_eof=True):
  if r.query_name in names and not r.is_secondary and not r.is_supplementary:
   rows.append(dict(name=r.query_name,mate=1 if r.is_read1 else 2,chrom=r.reference_name,start0=r.reference_start,end0=r.reference_end,reverse=r.is_reverse,mapq=r.mapping_quality,cigar=r.cigarstring,duplicate=r.is_duplicate,proper_pair=r.is_proper_pair))
(O/'supporting_read_pairs.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
