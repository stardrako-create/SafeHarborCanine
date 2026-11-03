from pathlib import Path
import json,pysam,datetime
A=Path(__file__).resolve().parent;O=A/'W11_ANCHOR_LOCALIZATION_20261002'
r=json.loads((A/'SV_PBGV000012/ROS_expanded_junction_20261001/review.json').read_text())
lengths=[]
for row in r['rows']:
 h=max(row['hits'],key=lambda h:h['AS']);assert h['cigar']=='150M'
 lengths.append({'read':row['read'],'left_aligned_bp':1000-h['start0'],'right_aligned_bp':h['end0']-1000})
maxright=max(x['right_aligned_bp'] for x in lengths)
with pysam.FastaFile(str(A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa')) as f:seq=f.fetch('junction_L4707939_R33613308',1000,1300).upper()
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results={}
for genome,path in [('ROS',A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'),('UU',Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa'))]:
 hits=json.loads((O/(genome+'_exact_hits.json')).read_text())['right_119'];ext=[]
 with pysam.FastaFile(str(path)) as f:
  for chrom,start,end,strand in hits:
   if strand=='+':target=f.fetch(chrom,start,min(start+300,f.get_reference_length(chrom))).upper()
   else:target=rc(f.fetch(chrom,max(0,end-300),end).upper())
   match=0
   for a,b in zip(seq,target):
    if a!=b:break
    match+=1
   ext.append({'chrom':chrom,'start0':start,'end0':end,'strand':strand,'exact_prefix_bp_capped300':match})
 results[genome]={'matches_at_max_observed_right_anchor':sum(x['exact_prefix_bp_capped300']>=maxright for x in ext),'loci':ext}
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'max_observed_right_aligned_bp':maxright,'read_anchor_lengths':lengths,'extension_from_all_exact119_hits':results,'limits':['Alignment lengths include possible mismatches; not all observed anchors are exact','Exact reference-prefix ambiguity does not exclude distinguishing variants in individual reads','Reference prefix over122bp is not established by these spanning reads','All longer exact hits must also contain the119bp prefix; enumerating these is exhaustive for exact-prefix comparison only']}
(O/'observed_anchor_extension.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'maxright':maxright,'counts':{k:v['matches_at_max_observed_right_anchor'] for k,v in results.items()},'ROS_loci':results['ROS']['loci']}))