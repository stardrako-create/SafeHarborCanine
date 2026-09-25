from pathlib import Path
import csv,json,subprocess,collections,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.1';O.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
 if not rows:return
 with p.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
windows=[r for r in read(P/'06_v1.22.0/frontier_T_evidence_v1220.tsv') if int(r['T_union_shared_fragment_starts'])>0];assert len(windows)==16
for i,r in enumerate(windows):r['window_id']=f'w{i:02d}';r['T_support_group']='consensus_and_union' if int(r['T_consensus_shared_fragment_starts']) else 'union_only'
write(O/'review_windows.tsv',windows)
ref=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';fa=pysam.FastaFile(str(ref));uu=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');bwa='/home/stardrako/miniforge3/envs/atac/bin/bwa';commands=[]
with (O/'windows.fa').open('w') as f:
 for r in windows:f.write(f'>{r["window_id"]}\n{fa.fetch(r["chrom"],int(r["start"]),int(r["end"]))}\n')
def align(index,query,out):
 cmd=[bwa,'mem','-a','-t','4',str(index),str(query)]
 with (O/(out+'.sam')).open('w') as f,(O/(out+'.log')).open('w') as e:p=subprocess.run(cmd,stdout=f,stderr=e)
 commands.append(dict(command=cmd,returncode=p.returncode));assert p.returncode==0
align(uu,O/'windows.fa','windows_uu')
svpath='/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz';sv=pysam.VariantFile(svpath);matches=[];variants=[]
with pysam.AlignmentFile(str(O/'windows_uu.sam'),'r') as sam:
 records=list(sam.fetch(until_eof=True));lengths=dict(zip(sam.references,sam.lengths))
 for r in windows:
  al=[a for a in records if a.query_name==r['window_id']];primary=[a for a in al if not a.is_secondary and not a.is_supplementary];assert len(primary)==1;a=primary[0];alt=len(al)>1 or a.has_tag('XA');nm=a.get_tag('NM') if a.has_tag('NM') else None
  aliases=[c for c,x in sv.header.contigs.items() if not a.is_unmapped and x.length==lengths[a.reference_name]]
  eligible=not a.is_unmapped and not alt and a.mapping_quality>=30 and a.query_alignment_length==1000 and nm is not None and nm<=10 and len(aliases)==1
  vname=aliases[0] if len(aliases)==1 else '';found=[]
  if eligible:
   for v in sv.fetch(vname,a.reference_start,a.reference_end):
    if v.start<a.reference_end and v.stop>a.reference_start:found.append(dict(window_id=r['window_id'],genes=r['genes'],chrom=vname,start=v.start,end=v.stop,id=v.id,SVTYPE=str(v.info.get('SVTYPE','')),AF=str(v.info.get('AF','')),filter=';'.join(v.filter)))
  matches.append(dict(window_id=r['window_id'],genes=r['genes'],ros_chrom=r['chrom'],ros_start=r['start'],ros_end=r['end'],uu_contig=a.reference_name,vcf_contig=vname,start=a.reference_start,end=a.reference_end,strand='-' if a.is_reverse else '+',mapq=a.mapping_quality,NM=nm,cigar=a.cigarstring,query_aligned_bp=a.query_alignment_length,alternatives=alt,eligible_for_variant_lookup=eligible,sv_records=len(found) if eligible else ''))
  variants.extend(found)
write(O/'dog10k_mapping.tsv',matches);write(O/'SV_records.tsv',variants)
tiles=[]
with (O/'tiles.fa').open('w') as f:
 for r in windows:
  for length in [100,150,250]:
   for s in range(int(r['start']),int(r['end'])-length+1,25):
    t=dict(tile_id=str(len(tiles)),window_id=r['window_id'],chrom=r['chrom'],start=s,end=s+length,length=length);tiles.append(t);f.write(f'>{t["tile_id"]}\n{fa.fetch(t["chrom"],s,s+length)}\n')
align(ref.with_suffix('.fa'),O/'tiles.fa','tiles_ros');als=collections.defaultdict(list)
with pysam.AlignmentFile(str(O/'tiles_ros.sam'),'r') as sam:
 for a in sam.fetch(until_eof=True):als[a.query_name].append(a)
tileout=[]
for t in tiles:
 al=als[t['tile_id']];pr=[a for a in al if not a.is_secondary and not a.is_supplementary];assert len(pr)==1;a=pr[0];correct=not a.is_unmapped and a.reference_name==t['chrom'] and a.reference_start==t['start'] and a.reference_end==t['end'] and a.get_tag('NM')==0
 other=any(not x.is_unmapped and (x.is_secondary or x.is_supplementary) for x in al) or a.has_tag('XA')
 tileout.append(dict(**t,mapq=a.mapping_quality,correct_exact_self_alignment=correct,other_alignment_reported=other,passes_remap_proxy=correct and a.mapping_quality>=30 and not other))
write(O/'tile_mappability_proxy.tsv',tileout);stats=[]
for r in windows:
 for length in [100,150,250]:
  rr=[t for t in tileout if t['window_id']==r['window_id'] and t['length']==length];stats.append(dict(window_id=r['window_id'],tile_bp=length,tested=len(rr),passed=sum(t['passes_remap_proxy'] for t in rr),min_mapq=min(t['mapq'] for t in rr)))
write(O/'tile_mappability_summary.tsv',stats)
(O/'mapping_execution.json').write_text(json.dumps(dict(commands=commands,windows=len(windows),eligible=sum(r['eligible_for_variant_lookup'] for r in matches),SV_records=variants,tile_tests=len(tiles),tile_proxy_failures=sum(not t['passes_remap_proxy'] for t in tileout),limitations=['BWA reported alignments are a mappability proxy, not exhaustive uniqueness or guide off-target assessment','Dog10K alias requires unique exact reference length match','16 overlapping windows are not independent loci','Reference differences do not identify the genotype of the experimental donor']),indent=2));print(json.dumps(dict(mapping=matches,tile_failures=[t for t in stats if t['passed']<t['tested']]),indent=2))
