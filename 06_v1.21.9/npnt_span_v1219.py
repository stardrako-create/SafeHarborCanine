from pathlib import Path
import pysam,json,csv,numpy as np,time
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.9'
review=json.loads((P/'06_v1.21.6/NPNT_SV_carriers_review.json').read_text());carriers=[r['sample'] for r in review['carriers']];meta={r['sampleName']:r for r in csv.DictReader((P/'06_v1.21.6/dog10K-alignment-sample-table.2022-02-23-v8.txt.UPDATESRA.txt').open(),delimiter='\t')}
url='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/AutoAndXPAR.SNPs.vqsr99.vcf.gz';vf=pysam.VariantFile(url)
controls=[s for s in vf.header.samples if s.startswith(('PBGV','GBGV','BFDB')) and s not in carriers];selected=carriers+controls;vf.subset_samples(selected)
intervals=[dict(name=f'inside_{i}',chrom='chr32',start=int(a),end=int(a)+2000,inside=True) for i,a in enumerate(np.linspace(6000000,32000000,10).astype(int))]+[dict(name=f'outside_{a}',chrom='chr32',start=a,end=a+2000,inside=False) for a in [1000000,3000000,35000000,37000000]]
rows=[];records=[]
for region in intervals:
 for r in vf.fetch(region['chrom'],region['start'],region['end']):
  if list(r.filter)!=['PASS']:continue
  for s in selected:
   g=r.samples[s];gt=g.get('GT');dp=g.get('DP');gq=g.get('GQ');called=gt is not None and all(a is not None for a in gt);het=called and len(set(gt))>1
   rows.append(dict(interval=region['name'],inside=region['inside'],chrom=r.chrom,pos=r.pos,ref=r.ref,alt=','.join(r.alts),sample=s,carrier=s in carriers,GT='/'.join('.' if a is None else str(a) for a in (gt or ())),DP=dp,GQ=gq,heterozygous=het,high_quality_het=bool(het and dp is not None and dp>=8 and gq is not None and gq>=20)))
 print(region['name'],'done',flush=True)
 with (O/'NPNT_span_genotypes.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary=[]
for s in selected:
 a=[r for r in rows if r['sample']==s and r['inside']];b=[r for r in rows if r['sample']==s and not r['inside']];inside=[r['DP'] for r in a if r['DP'] is not None];outside=[r['DP'] for r in b if r['DP'] is not None];im=float(np.median(inside)) if inside else None;om=float(np.median(outside)) if outside else None
 summary.append(dict(sample=s,carrier=s in carriers,breed=meta.get(s,{}).get('Breed/Type'),inside_sites=len(a),inside_high_quality_hets=sum(r['high_quality_het'] for r in a),inside_intervals_with_high_quality_het=len({r['interval'] for r in a if r['high_quality_het']}),inside_median_DP=im,outside_median_DP=om,inside_outside_DP_ratio=im/om if im is not None and om else None))
(O/'NPNT_span_review.json').write_text(json.dumps(dict(source=url,intervals=intervals,carrier_samples=carriers,breed_control_samples=controls,summary=summary,limitations=['Variant-ascertained sparse sites; DP is not unbiased genome-wide copy number','Calls cannot exclude mosaicism, complex rearrangements or mapping artefacts','No raw-read or orthogonal SV validation','High-quality heterozygosity inside an asserted deletion challenges a simple constitutive one-copy loss at those positions']),indent=2));print(json.dumps([r for r in summary if r['carrier']],indent=2))
