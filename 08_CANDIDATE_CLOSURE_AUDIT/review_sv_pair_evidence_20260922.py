"""Follow-up to the anchor-only screen; actual mate records and separate evidence classes."""
from pathlib import Path
import collections,datetime,json,pysam
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked'
CHROM='NC_049253.1'; BOUNDS=[4707938,33613308]
def classify(left,right):
 if left is None or right is None:return 'missing_or_ambiguous_mate_record'
 for r in (left,right):
  if r['mapq']<30:return 'low_MAPQ'
  if r['unmapped'] or r['secondary'] or r['supplementary'] or r['qcfail'] or r['duplicate'] or not r['paired']:return 'ineligible_record'
 if left['mate']==right['mate'] or left['mate'] not in (1,2) or right['mate'] not in (1,2):return 'mate_identity_mismatch'
 if left['chrom']!=CHROM or right['chrom']!=CHROM or abs(left['start0']-BOUNDS[0])>10000 or abs(right['start0']-BOUNDS[1])>10000:return 'outside_boundary_neighborhood'
 if any(r['mate_unmapped'] for r in (left,right)) or left['mate_chrom']!=right['chrom'] or right['mate_chrom']!=left['chrom'] or left['mate_start0']!=right['start0'] or right['mate_start0']!=left['start0']:return 'nonreciprocal_mate_fields'
 if left['mate_reverse']!=right['reverse'] or right['mate_reverse']!=left['reverse']:return 'inconsistent_mate_strand_fields'
 return 'both_MAPQ30_inward' if not left['reverse'] and right['reverse'] else 'both_MAPQ30_other_orientation'
def main():
 old=json.loads((A/'receipt_and_boundary_link_review_20260922.json').read_text())
 names=set().union(*(set(r['mate_links_to_other_boundary']) for r in old['boundary_link_screen']))
 found=collections.defaultdict(lambda:[[],[]])
 with pysam.AlignmentFile(str(A/'selected_markdup.bam'),'rb') as bam:
  for side,pos in enumerate(BOUNDS):
   for r in bam.fetch(CHROM,pos-10000,pos+10000):
    if r.query_name not in names or r.is_secondary or r.is_supplementary:continue
    found[r.query_name][side].append(dict(mate=1 if r.is_read1 and not r.is_read2 else 2 if r.is_read2 and not r.is_read1 else 0,chrom=r.reference_name,start0=r.reference_start,mapq=r.mapping_quality,reverse=r.is_reverse,paired=r.is_paired,unmapped=r.is_unmapped,secondary=r.is_secondary,supplementary=r.is_supplementary,qcfail=r.is_qcfail,duplicate=r.is_duplicate,mate_unmapped=r.mate_is_unmapped,mate_chrom=r.next_reference_name,mate_start0=r.next_reference_start,mate_reverse=r.mate_is_reverse))
 rows=[]
 for name in sorted(names):
  sides=found[name]; l,r=[x[0] if len(x)==1 else None for x in sides]
  rows.append(dict(name=name,classification=classify(l,r),left_records=sides[0],right_records=sides[1]))
 sa=sorted(set().union(*(set(r['SA_links_to_other_boundary']) for r in old['boundary_link_screen'])))
 out=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source='receipt_and_boundary_link_review_20260922.json',anchor_only_screen_unique_names=len(names),pair_record_review=rows,pair_class_counts=dict(collections.Counter(r['classification'] for r in rows)),SA_anchor30_tag30_unique_names=sa,limitations=['Follow-up restricted to candidates of the original anchor-only screen','SA evidence retains tag-based MAPQ filtering; no supplementary-record or junction-orientation validation claimed','Inward orientation is compatible with a deletion, not confirmation','Pair and SA evidence are separate; no combined support count; read names may overlap','Repeat ambiguity and coverage findings remain applicable'])
 (A/'boundary_link_symmetric_pair_review_20260922.json').write_text(json.dumps(out,indent=2))
 print(json.dumps({k:v for k,v in out.items() if k!='pair_record_review'},indent=2))
if __name__=='__main__':main()
