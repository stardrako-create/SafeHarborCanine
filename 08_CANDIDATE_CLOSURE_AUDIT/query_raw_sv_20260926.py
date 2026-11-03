from pathlib import Path
import pysam,json,urllib.request,datetime
A=Path(__file__).resolve().parent;O=A/'SV_CALL_AUDIT_20260926'
u='https://kiddlabshare.med.umich.edu/dog10K/Manta-SV_2022-03-28/raw/SV-genotype-v2.merge.vcf.gz'
idx=O/'raw_remote.vcf.gz.tbi'
if not idx.exists():urllib.request.urlretrieve(u+'.tbi',idx)
v=pysam.VariantFile(u,index_filename=str(idx));v.subset_samples(['PBGV000010']);rows=[]
for r in v.fetch('chr32',4707937,4707941):
 if r.pos<4707938 or r.pos>4707940:continue
 rows.append(dict(id=r.id,pos=r.pos,end=r.stop,alts=r.alts,info=dict(r.info),filters=list(r.filter),sample=dict(r.samples['PBGV000010'])))
(O/'raw_region_review.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'url':u,'rows':rows},indent=2));print('Regional records:',len(rows));print(json.dumps([{'id':r['id'],'end':r['end'],'alts':r['alts']} for r in rows]))
