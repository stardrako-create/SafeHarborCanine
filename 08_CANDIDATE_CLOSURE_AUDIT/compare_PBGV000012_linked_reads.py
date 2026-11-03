from pathlib import Path
import pysam,json,subprocess,datetime
A=Path(__file__).resolve().parent;O=A/'SV_PBGV000012/ROS_structural_screen_20261001';r=json.loads((O/'review.json').read_text());names={x['name'] for x in r['reciprocal_mate_links_both_MAPQ30']};reads={}
with pysam.AlignmentFile(r['bam'],'rb') as f:
 for s,e in r['regions0'].values():
  for x in f.fetch('NC_051836.1',s,e):
   if x.query_name in names and not(x.is_secondary or x.is_supplementary):reads[(x.query_name,1 if x.is_read1 else 2)]=x.get_forward_sequence()
assert len(reads)==4
with (O/'linked_reads.fa').open('w') as f:
 for (n,m),s in sorted(reads.items()):f.write(f'>{n}_mate{m}\n{s}\n')
ref=A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa'
assert Path(str(ref)+'.bwt').exists()
with (O/'joint.sam').open('w') as out,(O/'joint.stderr').open('w') as err:subprocess.run(['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','2',str(ref),str(O/'linked_reads.fa')],stdout=out,stderr=err,check=True)
hits={}
with pysam.AlignmentFile(str(O/'joint.sam'),'r') as f:
 for x in f:
  if not x.is_unmapped:hits.setdefault(x.query_name,[]).append({'ref':x.reference_name,'start0':x.reference_start,'end0':x.reference_end,'AS':x.get_tag('AS'),'mapq':x.mapping_quality,'cigar':x.cigarstring,'blocks':x.get_blocks()})
(O/'joint_hits.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'hits':hits,'scope':'four selected reads; single-read competition only, not independent validation'},indent=2))
print('joint comparison completed')