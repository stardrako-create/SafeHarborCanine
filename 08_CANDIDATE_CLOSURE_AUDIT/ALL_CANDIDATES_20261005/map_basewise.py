from pathlib import Path
import json,pysam,collections
O=Path(__file__).parent
out=[]
with pysam.AlignmentFile(str(O/'UU_projection/full_reference.sam'),'r') as sam,pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa') as fa:
 for r in sam:
  if r.is_unmapped or r.is_secondary or r.is_supplementary:continue
  seq=r.query_sequence.upper();pairs=r.get_aligned_pairs();blocks=[];mismatches=[];inserted=[];deleted=[];current=[]
  for q,t in pairs:
   if q is None:deleted.append(t)
   elif t is None:inserted.append(999-q if r.is_reverse else q)
   if q is not None and t is not None:
    ref=fa.fetch(r.reference_name,t,t+1).upper();same=seq[q]==ref
    if not same:mismatches.append({'ROS_offset0':999-q if r.is_reverse else q,'UU_position0':t,'query_in_UU_orientation':seq[q],'UU_base':ref})
    if same:
     if current and (q!=current[-1][0]+1 or t!=current[-1][1]+1):blocks.append(current);current=[]
     current.append((q,t));continue
   if current:blocks.append(current);current=[]
  if current:blocks.append(current)
  assert len(mismatches)+len(inserted)+len(deleted)==r.get_tag('NM')
  b=[dict(ROS_offset_start0=1000-x[-1][0]-1 if r.is_reverse else x[0][0],ROS_offset_end0=1000-x[0][0] if r.is_reverse else x[-1][0]+1,UU_start0=x[0][1],UU_end0=x[-1][1]+1,length=len(x)) for x in blocks]
  assert sum(x['length'] for x in b)+len(mismatches)+len(inserted)==1000
  out.append(dict(id=r.query_name,UU_chrom=r.reference_name,strand='-' if r.is_reverse else '+',mapq=r.mapping_quality,cigar=r.cigarstring,exact_blocks=b,mismatches=mismatches,ROS_inserted_offsets=inserted,UU_deleted_positions=deleted))
assert len(out)==64
(O/'basewise_projection.json').write_text(json.dumps(dict(results=out,limitations=['Exact blocks only; no whole-window clearance for gapped/mismatched alignments','Indels crossing a block boundary require allele/haplotype reconciliation','Primary alignment is not proof of exhaustive uniqueness']),indent=2))
print('Verified 64 basewise maps, NM and query-length accounting')