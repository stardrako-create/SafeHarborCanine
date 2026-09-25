from pathlib import Path
import json,subprocess,pysam,datetime
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked'; O=A/'junction_review_20260922';O.mkdir(exist_ok=True)
ref='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa'; chrom='NC_049253.1'; boundaries=[4707938,33613308]
rows=[];seen=set()
with pysam.AlignmentFile(str(A/'selected_markdup.bam'),'rb') as f:
 for side,pos in enumerate(boundaries):
  for r in f.fetch(chrom,pos-1000,pos+1000):
   if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate or r.mapping_quality<30:continue
   for end,ci in [('left',r.cigartuples[0]),('right',r.cigartuples[-1])]:
    op,n=ci
    if op!=4 or n<20:continue
    key=(r.query_name,r.is_read1,end)
    if key in seen:continue
    seen.add(key); seq=r.query_sequence[:n] if end=='left' else r.query_sequence[-n:]
    rows.append(dict(id='clip'+str(len(rows)),name=r.query_name,mate=1 if r.is_read1 else 2,side=side,clip_end=end,clip_length=n,anchor_start0=r.reference_start,anchor_end0=r.reference_end,anchor_cigar=r.cigarstring,anchor_reverse=r.is_reverse,sequence=seq))
(O/'clips.json').write_text(json.dumps(rows,indent=2));(O/'clips.fa').write_text(''.join('>'+r['id']+'\n'+r['sequence']+'\n' for r in rows))
cmd=['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','4',ref,str(O/'clips.fa')]
(O/'command.json').write_text(json.dumps(cmd))
with (O/'clips_full_reference.sam').open('w') as out,(O/'bwa.stderr.log').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True)
lookup={r['id']:r for r in rows};hits=[]
with pysam.AlignmentFile(str(O/'clips_full_reference.sam'),'r') as f:
 for r in f:
  if r.is_unmapped:continue
  original=lookup[r.query_name];other=boundaries[1-original['side']]
  hits.append(dict(clip_id=r.query_name,chrom=r.reference_name,start0=r.reference_start,end0=r.reference_end,mapq=r.mapping_quality,cigar=r.cigarstring,reverse=r.is_reverse,secondary=r.is_secondary,supplementary=r.is_supplementary,NM=r.get_tag('NM'),AS=r.get_tag('AS'),aligned_query_bases=r.query_alignment_length,near_opposite_boundary=r.reference_name==chrom and abs(r.reference_start-other)<=1000))
(O/'clip_hits.json').write_text(json.dumps(hits,indent=2)); links=[h for h in hits if h['near_opposite_boundary']]
summary=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),clips=len(rows),mapped_clips=len({h['clip_id'] for h in hits}),opposite_boundary_hits=links,limitations=['BWA mem default seed/score thresholds; short clips can be unmapped despite homology','Clip-only mapping quality is not an SV likelihood','Anchors within +/-1kb; clip length >=20; MAPQ30 nonduplicate primary anchors','Clips are aligned in BAM reference-oriented sequence; junction orientation must be interpreted with anchor'])
(O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
