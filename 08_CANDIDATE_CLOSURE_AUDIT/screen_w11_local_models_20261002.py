from pathlib import Path
import json,pysam,subprocess,collections,datetime,hashlib
A=Path(__file__).resolve().parent;I=A/'SV_PBGV000012/ROS_expanded_junction_20261001';O=A/'W11_LOCAL_MODELS_20261002';O.mkdir(exist_ok=False)
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
J=A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa'
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'
with pysam.FastaFile(str(J)) as f:j=f.fetch('junction_L4707939_R33613308').upper()
models={};meta={}
with pysam.FastaFile(str(R)) as f:
 for length in [164,181,208,300]:
  for b in [35342119,35342133]:
   key=f'local_repeat{length}_return{b}';models[key]=j[:1000]+j[1000:1000+length]+rc(f.fetch('NC_051836.1',b-1000,b).upper());meta[key]={'second_boundary0':1000+length,'repeat_length_hypothesis':length,'ROS_return_boundary0':b}
 models['ROS_reference_reverse']=rc(f.fetch('NC_051836.1',35341119,35343119).upper())
 models['catalogue_junction']=j
with (O/'models.fa').open('w') as f:
 for k,v in models.items():f.write('>'+k+'\n'+v+'\n')
(O/'models.json').write_text(json.dumps(meta,indent=2))
bwa='/home/stardrako/miniforge3/envs/atac/bin/bwa'
with (O/'index.log').open('w') as e:subprocess.run([bwa,'index',str(O/'models.fa')],stdout=e,stderr=e,check=True)
with (O/'screen.sam').open('w') as f,(O/'bwa.log').open('w') as e:subprocess.run([bwa,'mem','-a','-t','4',str(O/'models.fa'),str(I/'reads.fastq')],stdout=f,stderr=e,check=True)
baseline={}
with pysam.AlignmentFile(str(I/'joint.sam'),'r') as f:
 for x in f:
  if not x.is_unmapped:baseline[x.query_name]=max(baseline.get(x.query_name,-999),x.get_tag('AS'))
rows=[];pairhits=collections.defaultdict(list)
pairs=json.loads((A/'W11_EXTENSION_PAIRS_20261002/review.json').read_text());names={x.rsplit('_mate',1)[0] for x in pairs['candidate_reads']}
with pysam.AlignmentFile(str(O/'screen.sam'),'r') as f:
 for x in f:
  if x.is_unmapped:continue
  h={'read':x.query_name,'model':x.reference_name,'start0':x.reference_start,'end0':x.reference_end,'reverse':x.is_reverse,'AS':x.get_tag('AS'),'cigar':x.cigarstring,'MAPQ_local_only':x.mapping_quality}
  if x.query_name.rsplit('_mate',1)[0] in names:pairhits[x.query_name].append(h)
  if x.reference_name in meta:
   b=meta[x.reference_name]['second_boundary0']
   if any(s<=b-20 and e>=b+20 for s,e in x.get_blocks()):
    h['previous_joint_best_AS']=baseline.get(x.query_name);h['margin_to_previous_joint']=h['AS']-baseline[x.query_name] if x.query_name in baseline else None;rows.append(h)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'second_boundary_spanning_alignments':len(rows),'distinct_second_boundary_reads':len({x['read'] for x in rows}),'positive_margin_reads':len({x['read'] for x in rows if x['margin_to_previous_joint'] is not None and x['margin_to_previous_joint']>0}),'rows':rows,'limits':['Hypothesis screen of8local models, not exhaustive assembly','Repeat sequence taken from catalogue model; length164/181/208/300 is hypothetical','Return boundaries test no duplication or14bp overlap hypothesis, not proven structure','Local-only alignment MAPQ is not genomic specificity','Score comparisons across separate alignments are diagnostic, not likelihoods','Same local reads reused; failure does not exclude other lengths/sequences']}
(O/'review.json').write_text(json.dumps(summary,indent=2));(O/'candidate_pair_hits.json').write_text(json.dumps(pairhits,indent=2));(O/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),O/'models.fa',I/'reads.fastq',O/'screen.sam']},indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ['rows','limits']}))