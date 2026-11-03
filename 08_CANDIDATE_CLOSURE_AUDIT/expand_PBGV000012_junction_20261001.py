from pathlib import Path
import pysam,json,subprocess,datetime,fcntl,collections
A=Path(__file__).resolve().parent;O=A/'SV_PBGV000012/ROS_expanded_junction_20261001';O.mkdir(exist_ok=True)
l=open(O/'worker.lock','w');fcntl.flock(l,fcntl.LOCK_EX|fcntl.LOCK_NB)
B=A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam'
R=A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa'
def state(s,**kw):(O/'status.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':s,**kw},indent=2))
try:
 state('selecting_local_names_and_retained_mates')
 names=set();regions=[(35322119,35362119),(6382614,6442614)]
 with pysam.AlignmentFile(str(B),'rb') as f:
  for start,end in regions:
   for r in f.fetch('NC_051836.1',start,end):
    if not(r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail):names.add(r.query_name)
 reads={};metadata={}
 with pysam.AlignmentFile(str(B),'rb') as f:
  for r in f:
   if r.query_name not in names or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:continue
   k=r.query_name+('_mate1' if r.is_read1 else '_mate2');seq=r.get_forward_sequence();q=r.get_forward_qualities()
   assert seq is not None and q is not None and k not in reads
   reads[k]=(seq,q);metadata[k]={'original_MAPQ':r.mapping_quality,'ref':r.reference_name,'start0':r.reference_start,'cigar':r.cigarstring}
 with (O/'reads.fastq').open('w') as f:
  for k,(seq,q) in sorted(reads.items()):f.write('@'+k+'\n'+seq+'\n+\n'+''.join(chr(v+33) for v in q)+'\n')
 (O/'selection.json').write_text(json.dumps({'regions0':regions,'names':len(names),'reads':len(reads),'metadata':metadata,'selection':'primary nonduplicate QCpass; no MAPQ cutoff; mates present in retained BAM only'},indent=2))
 state('aligning_ROS_plus_nine_junctions',reads=len(reads))
 with (O/'joint.sam').open('w') as out,(O/'bwa.log').open('w') as err:subprocess.run(['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','6',str(R),str(O/'reads.fastq')],stdout=out,stderr=err,check=True)
 state('reviewing_competition')
 genomic={};spans=collections.defaultdict(list)
 with pysam.AlignmentFile(str(O/'joint.sam'),'r') as f:
  for r in f:
   if r.is_unmapped:continue
   score=r.get_tag('AS')
   if not r.reference_name.startswith('junction_'):genomic[r.query_name]=max(score,genomic.get(r.query_name,-100000))
   elif any(s<=980 and e>=1020 for s,e in r.get_blocks()):spans[r.query_name].append({'ref':r.reference_name,'AS':score,'MAPQ':r.mapping_quality,'cigar':r.cigarstring,'start0':r.reference_start,'end0':r.reference_end})
 rows=[]
 for name,hits in spans.items():
  best=max(h['AS'] for h in hits);g=genomic.get(name);rows.append({'read':name,'original':metadata[name],'best_spanning_AS':best,'best_reported_genomic_AS':g,'margin':best-g if g is not None else None,'hits':hits})
 summary={'reads_tested':len(reads),'names_selected':len(names),'spanning_read_count':len(rows),'spanning_with_positive_AS_margin':sum(r['margin'] is not None and r['margin']>0 for r in rows),'rows':rows,'limits':['Local windows plus mates retained in selected chr32 BAM; not all raw reads','Single-read competition, not calibrated pair likelihood','BWA reported hits are heuristic; absent genomic score is not infinite margin','Synthetic junction offset1000 and contiguous20bp anchors; ROS projection remains approximate','Same reads reused; not independent confirmation']}
 (O/'review.json').write_text(json.dumps(summary,indent=2));state('completed',reads=len(reads),spanning=len(rows));print(json.dumps({k:v for k,v in summary.items() if k not in ['rows','limits']}))
except Exception as e:state('failed',error=str(e));raise