from pathlib import Path
import pysam,csv,json,collections
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.6'
windows=list(csv.DictReader((P/'06_v1.21.5/dog10k_local_mapping.tsv').open(),delimiter='\t'))
fa=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');result=[];checked=0
for name in ['AutoAndXPAR.SNPs.vqsr99.vcf.gz','AutoAndXPAR.nonSNPs.filter.vcf.gz']:
 vf=pysam.VariantFile(str(O/(name+'.local.vcf')));records=list(vf);remote=json.loads((O/(name+'.lookup.json')).read_text())
 for w in windows:
  rows=[r for r in records if r.chrom==w['vcf_contig'] and r.start<int(w['end']) and r.stop>int(w['start'])]
  assert len(rows)==next(s['records'] for s in remote['summary'] if s['genes']==w['genes'])
  for r in rows:
   assert fa.fetch(w['uu_contig'],r.start,r.start+len(r.ref)).upper()==r.ref.upper(),str(r)
   checked+=1
  passed=[r for r in rows if list(r.filter)==['PASS']]
  af=lambda r: r.info.get('AF',()) if isinstance(r.info.get('AF',()),tuple) else (r.info.get('AF'),)
  result.append(dict(callset=name,genes=w['genes'],all_records=len(rows),PASS_records=len(passed),PASS_AF_ge_001=sum(any(a is not None and a>=.01 for a in af(r)) for r in passed),max_PASS_AF=max((a for r in passed for a in af(r) if a is not None),default=None),callset_samples=len(vf.header.samples)))
with (O/'small_variant_frequencies.tsv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');wr.writeheader();wr.writerows(result)
(O/'variant_validation.json').write_text(json.dumps(dict(reference_alleles_checked=checked,all_reference_alleles_match=True,independent_htslib_remote_counts_match=True,results=result,limitations=['1987-sample release; do not equate to published 1929-sample subset or 76 project dogs','AF is catalogue allele frequency, not individual carrier probability','SNP/nonSNP release and SV release have different sample sets','Window overlaps on UU assembly; no insertion point or guide selected']),indent=2))
print(json.dumps(result,indent=2))
