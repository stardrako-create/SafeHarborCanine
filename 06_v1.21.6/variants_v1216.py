from pathlib import Path
import pysam,csv,json,time
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'06_v1.21.6';O.mkdir(exist_ok=True)
windows=list(csv.DictReader((P/'06_v1.21.5/dog10k_local_mapping.tsv').open(),delimiter='\t'))
base='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/'
for name in ['AutoAndXPAR.SNPs.vqsr99.vcf.gz','AutoAndXPAR.nonSNPs.filter.vcf.gz']:
 print('Opening',name,flush=True)
 vf=pysam.VariantFile(base+name)
 results=[];summary=[]
 for w in windows:
  c=w['vcf_contig'];a=int(w['start']);b=int(w['end'])
  assert c in vf.header.contigs
  n=0
  for r in vf.fetch(c,a,b):
   if r.start>=b or r.stop<=a:continue
   n+=1
   results.append(dict(genes=w['genes'],chrom=c,start=r.start,end=r.stop,ref=r.ref,alt=','.join(r.alts or []),filter=';'.join(r.filter),AF=json.dumps(r.info.get('AF')),AC=json.dumps(r.info.get('AC')),AN=r.info.get('AN'),qual=r.qual))
  summary.append(dict(genes=w['genes'],records=n))
  print(name,w['genes'],n,flush=True)
 if results:
  with (O/(name+'.local.tsv')).open('w') as f:
   wr=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');wr.writeheader();wr.writerows(results)
 (O/(name+'.lookup.json')).write_text(json.dumps(dict(url=base+name,samples=len(vf.header.samples),windows=windows,summary=summary,header=str(vf.header)),indent=2))
 vf.close()
