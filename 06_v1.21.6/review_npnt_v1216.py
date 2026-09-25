from pathlib import Path
import pysam,csv,json
O=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/06_v1.21.6')
meta={r['sampleName']:r for r in csv.DictReader((O/'dog10K-alignment-sample-table.2022-02-23-v8.txt.UPDATESRA.txt').open(),delimiter='\t')}
vf=pysam.VariantFile('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz')
sv=next(r for r in vf.fetch('chr32',13062402,13063402) if r.id=='chr32:4707939:DG')
carriers={s:dict(v) for s,v in sv.samples.items() if any(a and a>0 for a in (v.get('GT') or ()))}
snps=pysam.VariantFile(str(O/'AutoAndXPAR.SNPs.vqsr99.vcf.gz.local.vcf'));records=[r for r in snps if r.chrom=='chr32' and list(r.filter)==['PASS']]
out=[]
for s,gt in carriers.items():
 het=[];called=0
 for r in records:
  if s not in r.samples:continue
  v=r.samples[s];g=v.get('GT')
  if g and all(a is not None for a in g):
   called+=1
   if len(set(g))>1:het.append(dict(pos=r.pos,GT=g,GQ=v.get('GQ'),DP=v.get('DP'),AD=v.get('AD')))
 out.append(dict(sample=s,metadata=meta.get(s),SV_genotype=gt,local_PASS_SNPs_called=called,local_heterozygous_SNPs=het))
(O/'NPNT_SV_carriers_review.json').write_text(json.dumps(dict(id=sv.id,info=dict(sv.info),carriers=out,interpretation='Catalogue genotype audit only. Heterozygous SNP evidence within a claimed deletion may be inconsistent with simple hemizygosity but requires read-level validation; no automatic acceptance or rejection.'),indent=2))
for r in out:print(r['sample'],r['metadata'].get('Breed/Type') if r['metadata'] else None,r['SV_genotype'].get('FT'),r['local_PASS_SNPs_called'],len(r['local_heterozygous_SNPs']),flush=True)
