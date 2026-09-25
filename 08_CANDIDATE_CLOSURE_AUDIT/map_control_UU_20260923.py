from pathlib import Path
import csv,json,pysam,subprocess,collections,datetime
P=Path(__file__).parent.parent;O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922/UU_projection';O.mkdir(exist_ok=True)
rows=list(csv.DictReader((O.parent/'transcript_and_boundary_review.tsv').open(),delimiter='\t'));R=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';U='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa';rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with pysam.FastaFile(str(R)) as fa:seqs={r['id']:fa.fetch(r['chrom'],int(r['start0']),int(r['end0'])).upper() for r in rows}
assert all(len(s)==1000 and not(set(s)-set('ACGT')) for s in seqs.values());(O/'controls_ROS.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in seqs.items()))
cmd=['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','4',U,str(O/'controls_ROS.fa')];(O/'command.json').write_text(json.dumps(cmd))
with (O/'full_reference.sam').open('w') as out,(O/'bwa.stderr.log').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True)
hits=collections.defaultdict(list)
with pysam.AlignmentFile(str(O/'full_reference.sam'),'r') as f,pysam.FastaFile(U) as ref:
 for r in f:
  if r.is_unmapped:continue
  exact=r.cigarstring=='1000M' and r.get_tag('NM')==0
  verified=False
  if exact:
   sequence=ref.fetch(r.reference_name,r.reference_start,r.reference_end).upper();verified=seqs[r.query_name]==(rc(sequence) if r.is_reverse else sequence);assert verified
  hits[r.query_name].append(dict(chrom=r.reference_name,start0=r.reference_start,end0=r.reference_end,cigar=r.cigarstring,NM=r.get_tag('NM'),mapq=r.mapping_quality,strand='-' if r.is_reverse else '+',primary=not r.is_secondary and not r.is_supplementary,exact_sequence_verified=verified))
results=[]
for key in seqs:
 exact=[r for r in hits[key] if r['exact_sequence_verified']];primary=[r for r in exact if r['primary'] and r['mapq']>=30]
 results.append(dict(id=key,status='eligible_exact_window_for_variant_query' if len(exact)==1 and len(primary)==1 else 'unresolved_for_simple_allele_projection',hits=hits[key]))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),results=results,counts=dict(collections.Counter(r['status'] for r in results)),limitations=['BWA-reported placements, not exhaustive genomic equivalence; independent of donor genotype','Eligibility only for full-window exact-sequence projection; variants spanning verified interval edges need separate treatment','No variant absence or AF conclusions from mapping alone'])
(O/'projection_review.json').write_text(json.dumps(result,indent=2));print(json.dumps(result['counts']))
