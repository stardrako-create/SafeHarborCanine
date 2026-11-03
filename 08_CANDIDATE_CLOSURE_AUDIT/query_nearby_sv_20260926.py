from pathlib import Path
import pysam,json,datetime
A=Path(__file__).resolve().parent;O=A/'SV_CALL_AUDIT_20260926';url='https://kiddlabshare.med.umich.edu/dog10K/Manta-SV_2022-03-28/raw/SV-genotype-v2.merge.vcf.gz'
v=pysam.VariantFile(url,index_filename=str(O/'raw_remote.vcf.gz.tbi'));v.subset_samples(['PBGV000010']);rows={}
for center in [4707939,33613308]:
 for r in v.fetch('chr32',center-10000,center+10000):
  gt=r.samples['PBGV000010'].get('GT');
  if not gt or not any(x is not None and x>0 for x in gt):continue
  if min(abs(r.pos-center),abs(r.stop-center))>10000:continue
  rows[r.id]=dict(id=r.id,pos=r.pos,end=r.stop,alts=r.alts,info=dict(r.info),filters=list(r.filter),sample=dict(r.samples['PBGV000010']))
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':url,'query':'10kb around both w11 catalogue endpoints; nonreference genotypes in PBGV000010','rows':list(rows.values()),'limits':'Unphased catalogue models; alternative calls can be correlated or artifacts, not independent observations'}
(O/'nearby_nonreference_calls.json').write_text(json.dumps(result,indent=2));print(json.dumps([{'id':r['id'],'pos':r['pos'],'end':r['end'],'GT':r['sample']['GT'],'FT':r['sample'].get('FT')} for r in rows.values()]))
