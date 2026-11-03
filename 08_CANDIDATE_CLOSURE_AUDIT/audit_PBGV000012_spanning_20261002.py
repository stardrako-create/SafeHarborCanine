from pathlib import Path
import json,pysam,subprocess,hashlib,datetime
A=Path(__file__).resolve().parent
I=A/'SV_PBGV000012/ROS_expanded_junction_20261001'
O=A/'SV_PBGV000012/ROS_spanning_audit_20261002';O.mkdir(exist_ok=False)
r=json.loads((I/'review.json').read_text()); keys={x['read'] for x in r['rows']}; names={k.rsplit('_mate',1)[0] for k in keys}
seqs={}
with pysam.FastxFile(str(I/'reads.fastq')) as f:
 for x in f:
  if x.name in keys:seqs[x.name]=x.sequence
assert len(seqs)==12
with (O/'spanning.fa').open('w') as f:
 for k,v in seqs.items():f.write('>'+k+'\n'+v+'\n')
bam=A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam'; records=[]
with pysam.AlignmentFile(str(bam),'rb') as f:
 for x in f:
  if x.query_name in names:
   records.append(dict(name=x.query_name,mate=1 if x.is_read1 else 2,ref=x.reference_name,start0=x.reference_start,end0=x.reference_end,reverse=x.is_reverse,mapq=x.mapping_quality,cigar=x.cigarstring,flag=x.flag,duplicate=x.is_duplicate,secondary=x.is_secondary,supplementary=x.is_supplementary,AS=x.get_tag('AS') if x.has_tag('AS') else None,SA=x.get_tag('SA') if x.has_tag('SA') else None))
(O/'original_records.json').write_text(json.dumps(records,indent=2))
ref=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'
with (O/'genomic.sam').open('w') as out,(O/'bwa.log').open('w') as err:subprocess.run(['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','4',str(ref),str(O/'spanning.fa')],stdout=out,stderr=err,check=True)
hits={k:[] for k in keys}
with pysam.AlignmentFile(str(O/'genomic.sam'),'r') as f:
 for x in f:
  if not x.is_unmapped:hits[x.query_name].append(dict(ref=x.reference_name,start0=x.reference_start,cigar=x.cigarstring,AS=x.get_tag('AS'),mapq=x.mapping_quality))
rows=[]
for x in r['rows']:
 k=x['read'];h=hits[k];g=max((t['AS'] for t in h),default=None)
 rows.append(dict(read=k,best_junction_AS=x['best_spanning_AS'],genomic_only_best_AS=g,margin=x['best_spanning_AS']-g if g is not None else None,genomic_hits=h))
signatures={}
for name in names:
 pp=[x for x in records if x['name']==name and not x['secondary'] and not x['supplementary']]
 sig=tuple(sorted((x['mate'],x['ref'],x['start0'],x['end0'],x['reverse'],x['cigar']) for x in pp))
 signatures.setdefault(str(sig),[]).append(name)
review=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),spanning_reads=len(keys),distinct_names=len(names),unique_sequences=len(set(seqs.values())),pair_alignment_signature_groups=list(signatures.values()),rows=rows,limits=['Separate genomic-only alignment is diagnostic, not joint likelihood or independent data','Different read names or pair signatures do not prove distinct original molecules','Repeat ambiguity and structural size remain unresolved; normal broad depth does not support simple large heterozygous deletion'])
(O/'review.json').write_text(json.dumps(review,indent=2))
(O/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),I/'review.json',I/'reads.fastq',O/'spanning.fa',O/'genomic.sam']},indent=2))
print(json.dumps(review))