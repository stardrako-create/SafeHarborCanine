from pathlib import Path
import csv,json,hashlib,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT'
fa=pysam.FastaFile(str(P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'))
source={r['id']:r for r in csv.DictReader((O/'variant_projection_provenance.tsv').open(),delimiter='\t')}
rows=list(csv.DictReader((O/'normalized_variant_alleles_ROS.tsv').open(),delimiter='\t'))
for r in rows:
 old=source[r['id']];pos=int(r['start0']);op=int(old['ROS_start0']);ref=r['ref'].upper();alt=r['alt'].upper()
 lo=min(pos,op)-50;hi=max(pos+len(ref),op+len(old['ROS_ref']))+50
 seq=fa.fetch(r['chrom'],lo,hi).upper()
 assert seq[pos-lo:pos-lo+len(ref)]==ref
 nh=seq[:pos-lo]+alt+seq[pos-lo+len(ref):]
 alts=old['ROS_alt'].upper().split(',');freq=json.loads(old['AF'])
 assert len(alts)==len(freq)
 matching=[i for i,a in enumerate(alts) if seq[:op-lo]+a+seq[op-lo+len(old['ROS_ref']):]==nh]
 assert len(matching)==1,(r['id'],matching)
 idx=matching[0];af=freq[idx];assert af is not None and 0<=af<=1
 r.update(source_alt_index_1based=idx+1,source_UU_alt=old['UU_alt'].split(',')[idx],allele_AF=af,allele_AF_ge_0_01=af>=0.01,allele_identity_method='unique_local_haplotype_equivalence')
dest=O/'normalized_variant_alleles_with_AF_ROS.tsv'
with dest.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary=[]
for wid in sorted({r['window_id'] for r in rows}):
 subset=[r for r in rows if r['window_id']==wid];passed=[r for r in subset if r['filter']=='PASS']
 summary.append(dict(window=wid,total_alleles=len(subset),PASS_alleles=len(passed),PASS_alleles_AF_ge_1pct=sum(r['allele_AF_ge_0_01'] for r in passed),PASS_alleles_misclassified_if_record_max_used=sum((float(r['source_record_max_alt_AF'])>=0.01)!=r['allele_AF_ge_0_01'] for r in passed)))
result=dict(status='completed',all_alleles_uniquely_matched=True,summary=summary,corrected_rows=[{k:r[k] for k in ['id','alt','source_record_max_alt_AF','allele_AF']} for r in rows if abs(float(r['source_record_max_alt_AF'])-r['allele_AF'])>1e-12],output_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),limitations=['Catalogue allele frequencies, not donor genotypes','No allele deleted or filtered; original outputs preserved','Normalization and strand reversal preserve original ALT ordering only through verified haplotype matching'])
(O/'allele_frequency_review.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
