"""Prepare reference-only coverage windows; no inference of sample copy number."""
from pathlib import Path
import csv,json,hashlib,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT';fa=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa')
chrom='NC_049253.1';length=fa.get_reference_length(chrom)
svs=list(csv.DictReader((P/'06_v1.22.1/SV_records.tsv').open(),delimiter='\t'))
sv=[r for r in svs if r['window_id']=='w11' and r['id']=='chr32:4707939:DG'];assert len(sv)==1
start,end=int(sv[0]['start']),int(sv[0]['end']);assert 0<start<end<length
rows=[]
for a in range(0,length,100000):
 b=min(a+100000,length);seq=fa.fetch(chrom,a,b).upper();acgt=sum(seq.count(c) for c in 'ACGT')
 label='inside' if start<=a and b<=end else 'outside' if b<=start or a>=end else 'crosses_catalogue_boundary'
 rows.append(dict(chrom=chrom,start0=a,end0=b,bases=b-a,relation_to_catalogued_span=label,GC_fraction_canonical=(seq.count('G')+seq.count('C'))/acgt if acgt else None,canonical_fraction=acgt/(b-a),N_bases=seq.count('N')))
with (O/'SV_chr32_100kb_reference_windows.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
with (O/'SV_chr32_100kb_windows.bed').open('w') as f:
 for r in rows:f.write(f"{chrom}\t{r['start0']}\t{r['end0']}\t{r['relation_to_catalogued_span']}\n")
with (O/'SV_breakpoint_neighborhoods_UU.bed').open('w') as f:
 for name,pos in [('left',start),('right',end)]:f.write(f'{chrom}\t{max(0,pos-10000)}\t{min(length,pos+10000)}\t{name}_catalogue_boundary_10kb_each_side\n')
assert sum(r['bases'] for r in rows)==length
assert all(a['end0']==b['start0'] for a,b in zip(rows,rows[1:]))
result={'status':'reference_intervals_prepared_not_sample_analysis','assembly':'UU_Cfam_GSD_1.0','chrom':chrom,'vcf_chrom_alias':'chr32','chrom_length':length,'catalogue_span_0based':[start,end],'windows':len(rows),'window_counts_by_relation':{k:sum(r['relation_to_catalogued_span']==k for r in rows) for k in ['inside','outside','crosses_catalogue_boundary']},'all_bases_covered_exactly_once':True,'limitations':['GC and N are reference descriptors, not sample coverage','No mappability correction yet; repeat-rich windows need separate assessment','Outside windows are candidate comparators, not validated diploid controls','Catalogue boundary coordinates retained exactly; breakpoint sequence itself has not been confirmed','No sample copy number or SV decision made']}
(O/'SV_coverage_window_preparation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
