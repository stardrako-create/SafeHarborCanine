from pathlib import Path
import json,pysam,collections,hashlib
O=Path(__file__).parent
pending=json.loads((O/'pending_classification_20261006.json').read_text())['observations']
records={};sources=[]
for folder in ['variant_catalogue_query_drop_samples','variant_catalogue_query_nonexact']:
 for f in sorted((O/folder).glob('*__*.vcf')):
  sources.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
  with pysam.VariantFile(str(f)) as vf:
   for r in vf:
    key=(r.contig,r.pos,r.ref,tuple(r.alts or ()))
    rec=dict(chrom=r.contig,pos1=r.pos,ref=r.ref,alts=list(r.alts or ()),AF=list(r.info.get('AF',())),filters=list(r.filter))
    if key in records:assert records[key]==rec
    records[key]=rec
stars=[];swaps=[]
for x in pending:
 a=x['allele'];pos=a['pos1']-1
 if a['alt']=='*':
  candidates=[]
  for r in records.values():
   if r['chrom']!=a['chrom']:continue
   for i,alt in enumerate(r['alts']):
    ref=r['ref'];start=r['pos1']-1
    if not alt or set(alt.upper())-set('ACGT') or len(ref)<=len(alt):continue
    while ref and alt and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
    while ref and alt and ref[0]==alt[0]:ref=ref[1:];alt=alt[1:];start+=1
    if not alt and start<=pos<start+len(ref):candidates.append(dict(record=r,alt_index1=i+1,deleted_start0=start,deleted_end0=start+len(ref)))
  stars.append(dict(window=x['window'],star=a,candidate_deletions=candidates,status='positional_candidates_not_phased' if candidates else 'not_explained_by_downloaded_records'))
 elif x['reason']=='catalogue_ALT_matches_ROS_reference':
  matches=[r for r in records.values() if r['chrom']==a['chrom'] and r['pos1']==a['pos1'] and r['ref']==a['ref'] and a['alt'] in r['alts']]
  assert len(matches)==1
  r=matches[0];assert len(r['AF'])==len(r['alts'])
  total=sum(r['AF']);assert -1e-6<=total<=1.00001
  swaps.append(dict(window=x['window'],record=r,ROS_reference_ALT_index1=a['alt_index1'],UU_reference_frequency_derived=max(0,1-total),method='1-sum(all original ALT AF), includes symbolic alleles; rounded source frequencies',not_individual_genotype=True))
report=dict(sources=sources,unique_records=len(records),star_observations=stars,reference_swaps=swaps,limitations=['Positional overlap cannot establish that a deletion and star occur on the same haplotype','Downloaded windows may omit upstream deletion records','Derived REF frequency is population aggregate with rounding, not donor genotype'])
(O/'star_and_reference_context_20261006.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'unique_records':len(records),'stars':len(stars),'stars_with_positional_candidates':sum(bool(x['candidate_deletions']) for x in stars),'reference_swaps':len(swaps)}))
w='NC_051831.1_39708074_39709074'
print(json.dumps({'priority_stars':[x for x in stars if x['window']==w],'priority_swaps':[x for x in swaps if x['window']==w]}))