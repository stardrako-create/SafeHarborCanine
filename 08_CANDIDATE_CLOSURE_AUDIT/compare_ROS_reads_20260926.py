from pathlib import Path
import pysam,json,subprocess,hashlib,datetime
A=Path(__file__).resolve().parent;P=A/'SV_PBGV000010/alignment_chunked';O=P/'ROS_read_comparison_20260926';O.mkdir(exist_ok=False)
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';reads=P/'whole_genome_competition_20260926/reads.fa';B='/home/stardrako/miniforge3/envs/atac/bin/bwa'
cmd=[B,'mem','-a','-t','4',str(R),str(reads)];(O/'command.json').write_text(json.dumps(cmd))
with (O/'ROS.sam').open('w') as f,(O/'stderr.log').open('w') as e:subprocess.run(cmd,stdout=f,stderr=e,check=True)
hits={}
with pysam.AlignmentFile(str(O/'ROS.sam')) as f:
 for x in f:
  if x.is_unmapped:continue
  hits.setdefault(x.query_name,[]).append(dict(chrom=x.reference_name,start0=x.reference_start,end0=x.reference_end,cigar=x.cigarstring,NM=x.get_tag('NM'),AS=x.get_tag('AS'),mapq=x.mapping_quality,reverse=x.is_reverse,secondary=x.is_secondary,supplementary=x.is_supplementary))
previous=json.loads((P/'whole_genome_competition_20260926/review.json').read_text())['results']['single'];rows=[]
for x in previous:
 hs=hits.get(x['read'],[]);best=max((h['AS'] for h in hs),default=None)
 rows.append(dict(read=x['read'],UU_best_genomic_AS=x['best_genomic_AS'],junction_AS=x['best_spanning_AS'],ROS_best_AS=best,ROS_best_hits=[h for h in hs if h['AS']==best],separate_run_junction_minus_ROS=None if best is None or x['best_spanning_AS'] is None else x['best_spanning_AS']-best))
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reads_sha256':hashlib.sha256(reads.read_bytes()).hexdigest(),'rows':rows,'limits':['ROS-only index versus prior UU-plus-junction index: score comparison descriptive, not joint ROS competition','Selected reads unchanged; not independent additional evidence','BWA reported alternatives are heuristic','No genotype or biological safe harbor validation']}
(O/'review.json').write_text(json.dumps(r,indent=2));(O/'all_hits.json').write_text(json.dumps(hits,indent=2));print(json.dumps([x for x in rows if x['junction_AS'] is not None],indent=2))
