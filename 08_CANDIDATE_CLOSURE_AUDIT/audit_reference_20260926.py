from pathlib import Path
import pysam,json,re,subprocess,datetime
A=Path(__file__).resolve().parent;O=A/'REFERENCE_AUDIT_20260926';O.mkdir(exist_ok=False)
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';U=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');B='/home/stardrako/miniforge3/envs/atac/bin/bwa'
header=(A/'W01_FOCUSED_REVIEW_20260923/SV_extended_20260924/header_without_samples.txt').read_text()
contigs=dict((n,int(v)) for n,v in re.findall(r'##contig=<ID=([^,]+),length=(\d+)',header))
with pysam.FastaFile(str(U)) as u,pysam.FastaFile(str(R)) as r:
 ul=dict(zip(u.references,u.lengths));rl=dict(zip(r.references,r.lengths));seq={}
 # Chromosome accessions ascend in both assemblies for autosomes 1-38.
 checks=[dict(chrom='chr'+str(i),catalogue=contigs.get('chr'+str(i)),UU=ul.get('NC_'+str(49221+i).zfill(6)+'.1'),ROS=rl.get('NC_'+str(51804+i).zfill(6)+'.1')) for i in range(1,39)]
 for label,pos in [('left',4707939),('right',33613308)]:
  seq[label+'_up']=u.fetch('NC_049253.1',pos-10000,pos);seq[label+'_down']=u.fetch('NC_049253.1',pos,pos+10000)
 seq['w11_UU']=u.fetch('NC_049253.1',13078664,13079664)
 seq['junction_best']=u.fetch('NC_049253.1',4706939,4707939)+u.fetch('NC_049253.1',33613308,33614308)
 (O/'assembly_lengths.json').write_text(json.dumps(checks,indent=2))
(O/'queries.fa').write_text(''.join('>'+n+'\n'+s+'\n' for n,s in seq.items()))
with (O/'queries_ROS.sam').open('w') as f,(O/'bwa.log').open('w') as e:subprocess.run([B,'mem','-a','-t','4',str(R),str(O/'queries.fa')],stdout=f,stderr=e,check=True)
hits=[]
with pysam.AlignmentFile(str(O/'queries_ROS.sam')) as f:
 for x in f:
  if not x.is_unmapped:hits.append(dict(query=x.query_name,chrom=x.reference_name,start0=x.reference_start,end0=x.reference_end,reverse=x.is_reverse,cigar=x.cigarstring,query_start=x.query_alignment_start,query_end=x.query_alignment_end,AS=x.get_tag('AS'),NM=x.get_tag('NM'),mapq=x.mapping_quality,secondary=x.is_secondary,supplementary=x.is_supplementary))
report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'explicit_reference_header':[l for l in header.splitlines() if l.startswith(('##reference=','##assembly='))],'autosomes_UU_lengths_match':sum(x['catalogue']==x['UU'] for x in checks),'autosomes_ROS_lengths_match':sum(x['catalogue']==x['ROS'] for x in checks),'hits':hits,'limits':'Length agreement supports assembly identification but not sequence identity; BWA flanking alignments are not reciprocal assembly validation.'}
(O/'review.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='hits'}));print(json.dumps([h for h in hits if not h['secondary']],indent=2))
