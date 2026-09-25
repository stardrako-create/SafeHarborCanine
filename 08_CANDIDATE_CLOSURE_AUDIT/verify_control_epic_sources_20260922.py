from pathlib import Path
import gzip,json,hashlib,datetime
P=Path(__file__).parent.parent;D=P/'06_v1.21.5/external_sources';O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922';manifest=json.loads((D/'epic_sources.json').read_text());out=[]
for s in manifest['sources']:
 p=D/Path(s['local']).name;h=hashlib.sha256();total=0
 with gzip.open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b);total+=len(b)
 out.append(dict(file=p.name,bytes_decompressed=total,sha256=h.hexdigest(),matches_manifest=h.hexdigest()==s['sha256_decompressed']))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),commit=manifest['commit'],files=out,all_match=all(r['matches_manifest'] for r in out),scope='Input integrity only; full canFam3.1 annotations. Regional ROS exports cover old candidates only; new control windows require audited cross-assembly mapping')
(O/'EpiC_source_integrity.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(files=len(out),all_match=result['all_match'],total_decompressed=sum(r['bytes_decompressed'] for r in out)),indent=2))
