from pathlib import Path
import json,pysam,re,collections,datetime,traceback
A=Path(__file__).resolve().parent;O=A/'W01_FOCUSED_REVIEW_20260923/SV_extended_20260924';O.mkdir(exist_ok=False)
source=Path('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz');prior=json.loads((O.parent/'focused_checks.json').read_text())
queries=[{'id':q['id'],'chrom':q['vcf_contig'],'start0':q['start0'],'end0':q['end0']} for q in prior['SV_overlap_query']];queries.append({'id':'w11_positive_context','chrom':'chr32','start0':13078664,'end0':13079664})
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
status={'state':'running','started_utc':now(),'records_scanned':0,'hits':[],'SVTYPE_counts':{},'BND_alt_count':0,'CIPOS_records':0,'CIEND_records':0,'source_bytes':source.stat().st_size,'queries':queries}
def save():
 p=O/'status.tmp';p.write_text(json.dumps(status,indent=2));p.replace(O/'status.json')
def overlaps(c,s,e,q):return c==q['chrom'] and s<q['end0'] and q['start0']<e
save();types=collections.Counter()
try:
 with pysam.VariantFile(str(source),drop_samples=True) as vf:
  status['INFO_fields']=list(vf.header.info);status['header_samples_ignored']=len(vf.header.samples)
  (O/'header_without_samples.txt').write_text(str(vf.header))
  for r in vf:
   status['records_scanned']+=1;types[str(r.info.get('SVTYPE','unknown'))]+=1
   ci1=r.info.get('CIPOS') if 'CIPOS' in vf.header.info else None;ci2=r.info.get('CIEND') if 'CIEND' in vf.header.info else None
   status['CIPOS_records']+=ci1 is not None;status['CIEND_records']+=ci2 is not None
   bounds=[]
   if ci1 is not None and len(ci1)==2 and None not in ci1:bounds.append(('CIPOS',r.contig,r.start+min(ci1),r.start+max(ci1)+1))
   if ci2 is not None and len(ci2)==2 and None not in ci2:bounds.append(('CIEND',r.contig,r.stop-1+min(ci2),r.stop-1+max(ci2)+1))
   mates=[]
   for alt in r.alts or ():
    match=re.search(r'[\[\]]([^\[\]]+):(\d+)[\[\]]',alt)
    if match:mates.append((match[1],int(match[2])-1));status['BND_alt_count']+=1
   for q in queries:
    evidence=[]
    if overlaps(r.contig,r.start,r.stop,q):evidence.append('record_interval_overlap')
    for label,c,s,e in bounds:
     if overlaps(c,s,e,q):evidence.append(label+'_overlap')
    for c,p in mates:
     if overlaps(c,p,p+1,q):evidence.append('BND_mate_position_overlap')
    if evidence:status['hits'].append({'query':q['id'],'contig':r.contig,'pos1':r.pos,'end0':r.stop,'id':r.id,'SVTYPE':r.info.get('SVTYPE'),'alts':r.alts,'filters':list(r.filter),'evidence':evidence})
   if status['records_scanned']%10000==0:status['SVTYPE_counts']=dict(types);status['updated_utc']=now();save()
 status['state']='completed';status['limitations']=['CIPOS/CIEND missing from catalogue cannot be recovered or treated as precise breakpoints','BND ALT mate positions screened globally; no sample genotype validation','Reported absence applies to this catalogue and these query windows only','w11 is a positive context check, not a resolved structural interpretation']
except Exception as e:status['state']='failed';status['error']=str(e);(O/'error.txt').write_text(traceback.format_exc())
finally:status['SVTYPE_counts']=dict(types);status['updated_utc']=now();save();print(json.dumps(status,indent=2))
