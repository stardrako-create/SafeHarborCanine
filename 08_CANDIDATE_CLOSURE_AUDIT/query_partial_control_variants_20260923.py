from pathlib import Path
import json,datetime,hashlib,signal,traceback
import pysam
A=Path(__file__).resolve().parent
O=A/'CONTROL_LOCUS_REVIEW_20260922/variant_catalogue_partial_query_20260923'
O.mkdir(exist_ok=True)
def save(name,obj):
 p=O/name; t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2));t.replace(p)
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def alarm(*args):raise TimeoutError('Remote operation exceeded 120 seconds')
signal.signal(signal.SIGALRM,alarm)
projection=A/'CONTROL_LOCUS_REVIEW_20260922/gapped_projection_20260923/UU_review.json'
windows=[r for r in json.loads(projection.read_text())['results'] if r['status']=='partial_mapping_only_not_window_clearance' and r['mapq']>=30 and r['additional_reported_alignments']==0]
base='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/'
status={'started_utc':now(),'state':'running','projection_sha256':hashlib.sha256(projection.read_bytes()).hexdigest(),'completed':[],'errors':[]}
save('status.json',status)
try:
 with pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa') as ref:
  lengths=dict(zip(ref.references,ref.lengths))
  for source in ['AutoAndXPAR.SNPs.vqsr99.vcf.gz','AutoAndXPAR.nonSNPs.filter.vcf.gz']:
   signal.alarm(120)
   vf=pysam.VariantFile(base+source)
   signal.alarm(0)
   (O/(source+'.header.txt')).write_text(str(vf.header))
   for w in windows:
    key=w['id']+'__'+source.replace('.vcf.gz','')
    if (O/(key+'.json')).exists():raise RuntimeError('Output already exists; refusing overwrite: '+key)
    h={'chrom':w['uu_chrom'],'start0':w['uu_start0'],'end0':w['uu_end0'],'strand':w['strand'],'partial_only':True}
    aliases=[c for c in vf.header.contigs if vf.header.contigs[c].length==lengths[h['chrom']]]
    if len(aliases)!=1:raise ValueError('Contig length alias ambiguous/missing: '+str((h['chrom'],aliases)))
    alias=aliases[0];records=[];alleles=[]
    signal.alarm(120)
    for r in vf.fetch(alias,h['start0'],h['end0']):
     records.append(str(r))
     observed=ref.fetch(h['chrom'],r.start,r.start+len(r.ref)).upper()
     ok=observed==r.ref.upper()
     af=r.info.get('AF')
     if af is not None and not isinstance(af,(tuple,list)):af=(af,)
     if af is not None and len(af)!=len(r.alts or ()):
      raise ValueError('AF cardinality does not match ALT count')
     for index,alt in enumerate(r.alts or ()):
      alleles.append({'chrom':alias,'uu_chrom':h['chrom'],'pos1':r.pos,'id':r.id,'ref':r.ref,'alt':alt,'alt_index1':index+1,'AF':None if af is None else af[index],'filters':list(r.filter),'reference_verified':ok,'within_queried_span_not_exact_window':r.start>=h['start0'] and r.start+len(r.ref)<=h['end0'],'sequence_allele':bool(alt) and set(alt.upper())<=set('ACGT')})
    signal.alarm(0)
    raw=O/(key+'.vcf');raw.write_text(str(vf.header)+''.join(records))
    result={'utc':now(),'source':base+source,'window':w['id'],'mapping':h,'vcf_contig':alias,'alias_method':'unique matching full contig length in this VCF header','record_count':len(records),'allele_count':len(alleles),'alleles':alleles,'raw_vcf_sha256':hashlib.sha256(raw.read_bytes()).hexdigest()}
    save(key+'.json',result)
    if any(not a['reference_verified'] for a in alleles):raise ValueError('REF mismatch: '+key)
    status['completed'].append({'key':key,'records':len(records),'alleles':len(alleles)})
    status['updated_utc']=now();save('status.json',status)
   vf.close()
 status['state']='query_completed_projection_pending'
except Exception as e:
 signal.alarm(0);status['state']='failed';status['errors'].append(str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:
 status['updated_utc']=now();save('status.json',status)
 print(json.dumps(status),flush=True)

