from pathlib import Path
import pysam,json,hashlib,datetime,traceback
O=Path(__file__).resolve().parent/'ANO2_star_genotype_context_20261009';O.mkdir(exist_ok=False)
base='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/'
queries=[('AutoAndXPAR.SNPs.vqsr99.vcf.gz',{6742368,6742773}),('AutoAndXPAR.nonSNPs.filter.vcf.gz',{6742362,6742771,6742772})];records=[];sources=[]
for filename,positions in queries:
 url=base+filename
 with pysam.VariantFile(url) as vf:
  dst=O/filename.removesuffix('.gz')
  with pysam.VariantFile(str(dst),'w',header=vf.header) as out:
   for r in vf.fetch('chr27',6742350,6742790):
    if r.pos not in positions:continue
    out.write(r);samples={}
    for name,call in r.samples.items():
     gt=call.get('GT');samples[name]=dict(GT=list(gt) if gt is not None else None,phased=call.phased,PS=call.get('PS'))
    records.append(dict(file=filename,chrom=r.contig,pos1=r.pos,ref=r.ref,alts=list(r.alts or ()),filters=list(r.filter),samples=samples))
 sources.append(dict(url=url,path=str(dst),sha256=hashlib.sha256(dst.read_bytes()).hexdigest()))
results=[]
for spos,dpos in [(6742368,6742362),(6742773,6742771),(6742772,6742771)]:
 stars=[r for r in records if r['pos1']==spos and '*' in r['alts']];dels=[r for r in records if r['pos1']==dpos and any(len(a)<len(r['ref']) and a!='*' for a in r['alts'])];assert len(stars)==len(dels)==1
 star=stars[0];deletion=dels[0];si=star['alts'].index('*')+1;di=next(i+1 for i,a in enumerate(deletion['alts']) if len(a)<len(deletion['ref']) and a!='*');carriers=[]
 for name,call in star['samples'].items():
  if si not in (call['GT'] or []):continue
  other=deletion['samples'][name];same_phase=call['phased'] and other['phased'] and call['PS'] is not None and call['PS']==other['PS']
  carriers.append(dict(sample=name,star=call,deletion=other,deletion_present=di in (other['GT'] or []),shared_explicit_phase_set=same_phase,matching_allele_slots=same_phase and [i for i,a in enumerate(call['GT']) if a==si]==[i for i,a in enumerate(other['GT']) if a==di]))
 results.append(dict(star_pos1=spos,deletion_pos1=dpos,star_carriers=len(carriers),carriers_with_deletion=sum(x['deletion_present'] for x in carriers),explicit_phase_matches=sum(x['matching_allele_slots'] for x in carriers),carriers=carriers))
receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=sources,results=results,limitations=['Same catalogue and samples, not independent confirmation.','Co-occurrence of unphased genotypes does not establish cis haplotype.','No individual experimental donor is required or inferred.']);(O/'review.json').write_text(json.dumps(receipt,indent=2));print(json.dumps([{k:v for k,v in x.items() if k!='carriers'} for x in results]))