from pathlib import Path
import pysam,csv,json,subprocess,collections
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');OLD=P/'06_v1.22.1';O=P/'06_v1.22.2';O.mkdir(exist_ok=True)
ros=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';uu=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');rf=pysam.FastaFile(str(ros));uf=pysam.FastaFile(str(uu));bwa='/home/stardrako/miniforge3/envs/atac/bin/bwa'
maps={r['window_id']:r for r in csv.DictReader((OLD/'dog10k_mapping.tsv').open(),delimiter='\t')};m=maps['w00'];rc=m['ros_chrom'];rs=int(m['ros_start']);re=int(m['ros_end']);uc=m['uu_contig'];us=int(m['start']);ue=int(m['end'])
with (O/'UU_reciprocal_queries.fa').open('w') as f:
 f.write('>w00_reciprocal\n'+uf.fetch(uc,us,ue)+'\n>w00_reciprocal_flanks\n'+uf.fetch(uc,us-1000,ue+1000)+'\n')
(O/'ROS_flank_query.fa').write_text('>w00_forward_flanks\n'+rf.fetch(rc,rs-1000,re+1000)+'\n');commands=[]
for ref,query,name in [(ros.with_suffix('.fa'),O/'UU_reciprocal_queries.fa','reciprocal_ros'),(uu,O/'ROS_flank_query.fa','flanks_uu')]:
 cmd=[bwa,'mem','-a','-t','4',str(ref),str(query)]
 with (O/(name+'.sam')).open('w') as out,(O/(name+'.log')).open('w') as err:p=subprocess.run(cmd,stdout=out,stderr=err)
 commands.append(dict(command=cmd,returncode=p.returncode));assert p.returncode==0;print(name,'done',flush=True)
with pysam.AlignmentFile(str(OLD/'windows_uu.sam'),'r') as f:original=next(a for a in f if a.query_name=='w00' and not a.is_secondary and not a.is_supplementary)
assert not original.is_reverse
forward={(rs+q,r) for q,r in original.get_aligned_pairs() if q is not None and r is not None};mismatches=[dict(ros_pos0=q,uu_pos0=r,ros_base=rf.fetch(rc,q,q+1),uu_base=uf.fetch(uc,r,r+1)) for q,r in sorted(forward) if rf.fetch(rc,q,q+1).upper()!=uf.fetch(uc,r,r+1).upper()]
aligns=[];backpairs=set()
for name in ['reciprocal_ros','flanks_uu']:
 groups=collections.defaultdict(list)
 with pysam.AlignmentFile(str(O/(name+'.sam')),'r') as f:
  for a in f:groups[a.query_name].append(a)
 for q,records in groups.items():
  primary=[a for a in records if not a.is_secondary and not a.is_supplementary];assert len(primary)==1;a=primary[0];other=len(records)>1 or a.has_tag('XA')
  aligns.append(dict(query=q,reference=a.reference_name,start=a.reference_start,end=a.reference_end,mapq=a.mapping_quality,NM=a.get_tag('NM'),cigar=a.cigarstring,query_aligned_bp=a.query_alignment_length,query_length=a.query_length,reverse=a.is_reverse,other_alignment_reported=other))
  if q=='w00_reciprocal':backpairs={(r,us+x) for x,r in a.get_aligned_pairs() if x is not None and r is not None};assert a.reference_name==rc and a.reference_start==rs and a.reference_end==re and not a.is_reverse
result=dict(original=m,forward_matched_or_mismatched_pair_bp=len(forward),ROS_unpaired_bp=1000-len(forward),UU_unpaired_bp=ue-us-len(forward),base_substitutions=mismatches,reciprocal_alignments=aligns,reciprocal_pair_intersection=len(forward&backpairs),forward_only_pair_count=len(forward-backpairs),reciprocal_only_pair_count=len(backpairs-forward),commands=commands,limitations=['Reciprocal alignments support reference correspondence, not individual donor genotype','BWA alignment uniqueness is heuristic; sequence within indels is not one-to-one','Supplementary regional variant lookup can use verified UU span while distinguishing unprojectable alleles from exact ROS coordinates','Original NM<=10 automated gate remains recorded; no biological veto inferred'])
assert len(mismatches)+result['ROS_unpaired_bp']+result['UU_unpaired_bp']==int(m['NM'])
(O/'LOC_mapping_resolution.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
# Preserve every reported alternative of failed NPNT tiles, including exact genomic positions.
tiles={r['tile_id']:r for r in csv.DictReader((OLD/'tile_mappability_proxy.tsv').open(),delimiter='\t') if r['passes_remap_proxy']=='False'};records=[]
with pysam.AlignmentFile(str(OLD/'tiles_ros.sam'),'r') as f:
 for a in f:
  if a.query_name in tiles:records.append(dict(tile_id=a.query_name,source_start=tiles[a.query_name]['start'],source_end=tiles[a.query_name]['end'],source_window=tiles[a.query_name]['window_id'],reference=a.reference_name,start=a.reference_start,end=a.reference_end,mapq=a.mapping_quality,cigar=a.cigarstring,NM=a.get_tag('NM') if a.has_tag('NM') else None,secondary=a.is_secondary,supplementary=a.is_supplementary,reverse=a.is_reverse,XA=a.get_tag('XA') if a.has_tag('XA') else ''))
with (O/'NPNT_failed_tile_alignments.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
assert len(tiles)==8
(O/'NPNT_mapping_review.json').write_text(json.dumps(dict(failed_tiles=8,source_footprint_start=min(int(t['start']) for t in tiles.values()),source_footprint_end=max(int(t['end']) for t in tiles.values()),alignment_records=records,limitations=['Footprint is union extent of tested failed segments, not exact minimal repeat boundary','Avoiding this footprint has not been shown to preserve T accessibility or editing specificity','No new window selected by clipping problematic bases']),indent=2));print('NPNT alternatives recorded',len(records))
