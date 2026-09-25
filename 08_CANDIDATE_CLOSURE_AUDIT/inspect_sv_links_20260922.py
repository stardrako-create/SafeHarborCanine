exec(open(__file__.replace('inspect_sv_links_20260922.py','review_sv_completion_20260922.py')).read().split('receipts=[]')[0])
review=json.loads((A/'receipt_and_boundary_link_review_20260922.json').read_text()); names=set(review['boundary_link_screen'][0]['mate_links_to_other_boundary']); rows=[]
with pysam.AlignmentFile(str(A/'selected_markdup.bam'),'rb') as f:
 for pos in [4707938,33613308]:
  for r in f.fetch('NC_049253.1',pos-10000,pos+10000):
   if r.query_name in names: rows.append(dict(name=r.query_name,flag=r.flag,start0=r.reference_start,end0=r.reference_end,mapq=r.mapping_quality,cigar=r.cigarstring,mate_start0=r.next_reference_start,template_length=r.template_length,SA=r.get_tag('SA') if r.has_tag('SA') else None))
(A/'boundary_link_alignments_20260922.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
