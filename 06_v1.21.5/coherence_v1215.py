from pathlib import Path
import csv,json,shutil,subprocess,os
import pyBigWig,numpy as np,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.5';BIN='/home/stardrako/miniforge3/envs/atac/bin/'
os.environ['PATH']=BIN+os.pathsep+os.environ['PATH']
bw=pyBigWig.open(str(P/'06_v1.20.0/ATAC/full76/peak_frequency.bw'))
bed=O/'consensus_peaks_76_min39_gap75.bed';n=0
with bed.open('w') as f:
 for ch in bw.chroms():
  blocks=[]
  for s,e,v in bw.intervals(ch) or []:
   if v<39:continue
   if blocks and s-blocks[-1][1]<=75:blocks[-1][1]=max(e,blocks[-1][1])
   else:blocks.append([s,e])
  for s,e in blocks:f.write(f'{ch}\t{s}\t{e}\n');n+=1
print('Consistent 76-dog consensus:',n,flush=True)
code=O/'code';code.mkdir(exist_ok=True)
for name in ['bw_utils.py','score_ship_candidates_v2.py']:shutil.copy2(P/'scripts'/name,code/name)
m=json.loads((P/'06_v1.21.1/manifest.json').read_text());cmd=m['command'];cmd[1]=str(code/'score_ship_candidates_v2.py')
for arg,value in {'--atac-peaks-bed':bed,'--out-scored':O/'candidates_scored_coherent76.tsv','--out-passing-bed':O/'candidates_passing_coherent76.bed'}.items():cmd[cmd.index(arg)+1]=str(value)
(O/'scoring_command.json').write_text(json.dumps(cmd,indent=2))
if not (O/'candidates_scored_coherent76.tsv').exists():
 with (O/'scoring.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
print('Coherent scorer completed',flush=True)
with (O/'cohort_and_support.tsv').open() as f:windows=[r for r in csv.DictReader(f,delimiter='\t') if not r['genes'].endswith('_region')]
ref=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';fa=pysam.FastaFile(str(ref));tiles=[]
with (O/'local_remap_tiles.fa').open('w') as f:
 for r in windows:
  for length in [100,150,250]:
   for s in range(int(r['start']),int(r['end'])-length+1,25):
    ident=str(len(tiles));tiles.append({'id':ident,'genes':r['genes'],'chrom':r['chrom'],'start':s,'end':s+length,'length':length});f.write(f'>{ident}\n{fa.fetch(r["chrom"],s,s+length)}\n')
index=ref.with_suffix('.fa')
with (O/'local_remap_tiles.sam').open('w') as out,(O/'local_remap_tiles.log').open('w') as err:subprocess.run(['bwa','mem','-a','-t','4',str(index),str(O/'local_remap_tiles.fa')],stdout=out,stderr=err,check=True)
aligns={str(i):[] for i in range(len(tiles))}
with pysam.AlignmentFile(str(O/'local_remap_tiles.sam'),'r') as sam:
 for a in sam.fetch(until_eof=True):aligns[a.query_name].append(a)
results=[]
for t in tiles:
 records=aligns[t['id']];primary=[a for a in records if not a.is_secondary and not a.is_supplementary];assert len(primary)==1
 a=primary[0];correct=not a.is_unmapped and a.reference_name==t['chrom'] and a.reference_start==t['start'] and a.reference_end==t['end'] and a.get_tag('NM')==0
 other=any(not x.is_unmapped and (x.is_secondary or x.is_supplementary) for x in records) or a.has_tag('XA')
 results.append({**t,'mapq':a.mapping_quality,'correct_exact_self_alignment':correct,'other_alignment_reported':other,'passes_remap_proxy':correct and a.mapping_quality>=30 and not other})
with (O/'local_remap_tiles.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
summary=[]
for gene in sorted({r['genes'] for r in results}):
 for length in [100,150,250]:
  rs=[r for r in results if r['genes']==gene and r['length']==length];summary.append({'genes':gene,'tile_bp':length,'tested':len(rs),'passes_remap_proxy':sum(r['passes_remap_proxy'] for r in rs),'min_mapq':min(r['mapq'] for r in rs)})
(O/'local_remap_summary.json').write_text(json.dumps(summary,indent=2))
shutil.copy2(__file__,O/'coherence_v1215.py')
print(json.dumps(summary,indent=2),flush=True)
