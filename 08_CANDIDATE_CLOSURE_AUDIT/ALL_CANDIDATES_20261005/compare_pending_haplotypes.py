from pathlib import Path
import json,pysam,collections
O=Path(__file__).parent
pending=json.loads((O/'pending_classification_20261006.json').read_text())['observations'];maps={}
with pysam.AlignmentFile(str(O/'UU_projection/full_reference.sam'),'r') as sam:
 for r in sam:
  if not r.is_unmapped and not r.is_secondary and not r.is_supplementary:maps[r.query_name]=(r.is_reverse,{t:q for q,t in r.get_aligned_pairs() if q is not None and t is not None})
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
results=[]
with pysam.FastaFile(str(O/'UU_projection/controls_ROS.fa')) as ros,pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa') as uu:
 for x in pending:
  if x['reason'] not in ['allele_crosses_indel_or_exact_block_boundary','allele_overlaps_reference_mismatch','crosses_alignment_envelope']:continue
  a=x['allele'];s=a['pos1']-1;e=s+len(a['ref']);rev,m=maps[x['window']];lo=s-50;hi=e+50
  result=dict(window=x['window'],allele=a,original_reason=x['reason'],status='flank_endpoints_unavailable')
  if lo in m and hi-1 in m:
   seq=ros.fetch(x['window']).upper();seq=rc(seq) if rev else seq
   qlo=m[lo];qhi=m[hi-1]+1
   assert qlo<qhi
   u=uu.fetch(a['uu_chrom'],lo,hi).upper();assert u[50:50+len(a['ref'])]==a['ref'].upper()
   hap=u[:50]+a['alt'].upper()+u[50+len(a['ref']):];rseq=seq[qlo:qhi]
   flank_ok=rseq[:50]==u[:50] and rseq[-50:]==u[-50:]
   result.update(status='catalogue_ALT_haplotype_equals_ROS' if flank_ok and hap==rseq else 'local_haplotypes_differ',flanks_50bp_exact=flank_ok,UU_interval=[lo,hi],ROS_oriented_interval=[qlo,qhi],UU_reference_sequence=u,UU_alternate_sequence=hap,ROS_sequence=rseq)
  results.append(result)
report=dict(results=results,counts=dict(collections.Counter(x['status'] for x in results)),limitations=['Exact local sequence equivalence only, not donor genotype','Other alleles/counts must be rebased separately','A mismatch is unresolved, not a biological exclusion','50bp flanks constrained to existing alignment; no new uniqueness assertion'])
(O/'indel_haplotype_comparison_20261008.json').write_text(json.dumps(report,indent=2));print(json.dumps(report['counts']))