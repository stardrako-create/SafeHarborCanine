from pathlib import Path
import csv,json,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.1'
windows=list(csv.DictReader((O/'dog10k_mapping.tsv').open(),delimiter='\t'));fa=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');base='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/';summary=[];allrows=[];provenance=[]
for name in ['AutoAndXPAR.SNPs.vqsr99.vcf.gz','AutoAndXPAR.nonSNPs.filter.vcf.gz']:
 vf=pysam.VariantFile(base+name);regions=[];cached={};count=0
 for c in sorted({w['vcf_contig'] for w in windows if w['eligible_for_variant_lookup']=='True'}):
  ws=[w for w in windows if w['vcf_contig']==c and w['eligible_for_variant_lookup']=='True'];a=min(int(w['start']) for w in ws);b=max(int(w['end']) for w in ws);regions.append(dict(chrom=c,start=a,end=b))
  cached[c]=list(vf.fetch(c,a,b));print(name,c,len(cached[c]),flush=True)
 with (O/(name+'.regional.vcf')).open('w') as f:
  f.write(str(vf.header))
  for records in cached.values():
   for r in records:f.write(str(r));count+=1
 for w in windows:
  eligible=w['eligible_for_variant_lookup']=='True';rr=[r for r in cached.get(w['vcf_contig'],[]) if eligible and r.start<int(w['end']) and r.stop>int(w['start'])];passed=[r for r in rr if list(r.filter)==['PASS']];common=[]
  for r in rr:
   observed=fa.fetch(w['uu_contig'],r.start,r.start+len(r.ref)).upper();assert observed==r.ref.upper(),(w['window_id'],r.pos,'REF mismatch')
   af=r.info.get('AF');afs=af if isinstance(af,tuple) else (af,);maxaf=max([x for x in afs if x is not None],default=None)
   if list(r.filter)==['PASS'] and maxaf is not None and maxaf>=.01:common.append(r)
   allrows.append(dict(window_id=w['window_id'],callset=name,chrom=r.chrom,start=r.start,end=r.stop,ref=r.ref,alt=','.join(r.alts or []),filter=';'.join(r.filter),AF=json.dumps(af),max_alt_AF=maxaf,AC=json.dumps(r.info.get('AC')),AN=r.info.get('AN'),REF_matches=True))
  summary.append(dict(window_id=w['window_id'],callset=name,eligible=eligible,records=len(rr) if eligible else '',PASS=len(passed) if eligible else '',PASS_max_alt_AF_ge_0_01=len(common) if eligible else ''))
 provenance.append(dict(url=base+name,samples=len(vf.header.samples),regions=regions,regional_records=count));vf.close()
for filename,rows in [('small_variant_records.tsv',allrows),('small_variant_summary.tsv',summary)]:
 with (O/filename).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(O/'variant_query_provenance.json').write_text(json.dumps(dict(sources=provenance,all_window_record_REF_checks_passed=True,limitations=['Window record counts overlap and are not unique population variants','AF is alternate allele frequency; no absence guarantee or experimental donor genotype','Ineligible mapping yields missing result, never zero','Regional VCFs preserve source genotypes and FILTER fields; no new variant calls']),indent=2));print(json.dumps(summary,indent=2))
