from pathlib import Path
import pysam,json,subprocess,datetime
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked';O=A/'junction_competition_20260922';O.mkdir(exist_ok=True)
R='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa';B='/home/stardrako/miniforge3/envs/atac/bin/bwa';C='NC_049253.1'
clips=json.loads((A/'junction_review_20260922/clips.json').read_text());names={r['name'] for r in clips};reads={}
with pysam.AlignmentFile(str(A/'selected_markdup.bam'),'rb') as f:
 for r in f.fetch(until_eof=True):
  if r.query_name in names and not r.is_secondary and not r.is_supplementary:
   key=r.query_name+('/1' if r.is_read1 else '/2');reads[key]=r.get_forward_sequence()
(O/'reads.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in sorted(reads.items())))
seqs={};meta={}
with pysam.FastaFile(R) as ref:
 for dl in [-1,0,1]:
  for dr in [-1,0,1]:
   l=4707938+dl;r=33613308+dr;key=f'junction_L{l}_R{r}'
   seqs[key]=ref.fetch(C,l-1000,l)+ref.fetch(C,r,r+1000);meta[key]=dict(left_end0=l,right_start0=r,junction_offset=1000)
 for label,p in [('left',4707938),('right',33613308)]:seqs['reference_'+label]=ref.fetch(C,p-1000,p+1000)
(O/'targets.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in seqs.items()));(O/'target_coordinates.json').write_text(json.dumps(meta,indent=2))
commands=[[B,'index',str(O/'targets.fa')],[B,'mem','-a','-t','4',R,str(O/'reads.fa')],[B,'mem','-a','-t','4',str(O/'targets.fa'),str(O/'reads.fa')]]
(O/'commands.json').write_text(json.dumps(commands,indent=2))
for cmd,output in zip(commands,['index.stdout','full_reference.sam','junction_targets.sam']):
 with (O/output).open('w') as out,(O/(output+'.stderr')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True)
full={};detail=[]
for source in ['full_reference','junction_targets']:
 with pysam.AlignmentFile(str(O/(source+'.sam')),'r') as f:
  for r in f:
   if r.is_unmapped:continue
   rec=dict(read=r.query_name,target=r.reference_name,start0=r.reference_start,end0=r.reference_end,AS=r.get_tag('AS'),NM=r.get_tag('NM'),cigar=r.cigarstring,mapq=r.mapping_quality,reverse=r.is_reverse,source=source)
   if source=='full_reference':full[r.query_name]=max(full.get(r.query_name,-999),rec['AS'])
   else:
    rec['left_aligned_bases']=sum(max(0,min(e,1000)-s) for s,e in r.get_blocks() if s<1000)
    rec['right_aligned_bases']=sum(max(0,e-max(s,1000)) for s,e in r.get_blocks() if e>1000)
    rec['junction_spanning']=r.reference_name.startswith('junction_') and min(rec['left_aligned_bases'],rec['right_aligned_bases'])>=20
   detail.append(rec)
spans=[dict(r,best_full_reference_AS=full.get(r['read']),score_gain_over_full_reference=r['AS']-full.get(r['read'],-999)) for r in detail if r.get('junction_spanning')]
(O/'all_alignments.json').write_text(json.dumps(detail,indent=2));(O/'junction_spanning_alignments.json').write_text(json.dumps(spans,indent=2))
best={}
for r in spans:
 if r['read'] not in best or r['AS']>best[r['read']]['AS']:best[r['read']]=r
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),selected_names=len(names),reads=len(reads),junctions_tested=9,best_spanning_per_read=list(best.values()),limitations=['Selected boundary-clipped reads and retained mates only; not exhaustive novel-SV discovery','Independent single-read BWA alignments, not paired likelihood or competitive joint index','BWA heuristic output and local scores are descriptive; synthetic-target MAPQ not comparable to genome MAPQ','A score advantage for a repetitive junction does not establish the genomic origin or variant genotype'])
(O/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
