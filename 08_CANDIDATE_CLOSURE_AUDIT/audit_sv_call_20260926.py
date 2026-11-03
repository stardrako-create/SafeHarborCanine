from pathlib import Path
import pysam,json,datetime,collections
A=Path(__file__).resolve().parent;O=A/'SV_CALL_AUDIT_20260926';O.mkdir(exist_ok=False)
v=pysam.VariantFile('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz')
samples=[x for x in v.header.samples if 'PBGV000010' in x];assert len(samples)==1,samples
rows=[]
for r in v.fetch('chr32',4707938,4707940):
 if r.pos!=4707939 or r.stop<33000000:continue
 counts=collections.Counter('/'.join('.' if a is None else str(a) for a in s.get('GT',())) for s in r.samples.values())
 rows.append(dict(id=r.id,pos1=r.pos,end1=r.stop,ref=r.ref,alts=r.alts,qual=r.qual,filters=list(r.filter),info=dict(r.info),sample=samples[0],sample_format=dict(r.samples[samples[0]]),genotype_counts=dict(counts),format_definitions={k:v.header.formats[k].description for k in r.format.keys()}))
assert len(rows)==1,len(rows)
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'records':rows,'limitations':['Catalogue genotype is a call, not independent confirmation','Population genotype totals are not allele frequency recalculation or independent validation','This sample is the public WGS individual, not an experimental donor']}
(O/'review.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
