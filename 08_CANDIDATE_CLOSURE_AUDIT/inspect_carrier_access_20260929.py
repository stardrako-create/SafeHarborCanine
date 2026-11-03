from pathlib import Path
import pysam,json,datetime,urllib.request,re,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_CARRIER_ACCESS_20260929';O.mkdir(exist_ok=False)
v=pysam.VariantFile('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz');rows=[]
for r in v.fetch('chr32',4707938,4707940):
 if r.id!='chr32:4707939:DG':continue
 for name,s in r.samples.items():
  if any(x is not None and x>0 for x in s.get('GT',())):rows.append({'sample':name,'format':dict(s)})
assert len(rows)==4
u='https://kiddlabshare.med.umich.edu/dog10K/cram-share/'
try:
 text=urllib.request.urlopen(u,timeout=25).read().decode();(O/'cram_directory.html').write_text(text);links=re.findall(r'href="([^"]+)"',text);matches={r['sample']:[x for x in links if r['sample'] in x] for r in rows};error=None
except Exception as e:matches={};error=type(e).__name__
report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'carriers':rows,'directory_url':u,'matching_links':matches,'access_error':error,'limits':['Catalogue carriers are hypotheses, not confirmed SV genotypes','FT PASS and failed calls must remain distinct','Availability of CRAM links does not confirm accessible indexed regional reads']}
(O/'review.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
