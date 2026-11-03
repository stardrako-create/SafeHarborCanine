from pathlib import Path
import csv,json,hashlib
N=Path(__file__).resolve().parent;P=N.parent.parent;src=P/'06_v1.21.8/candidates_integrated_review_v1218.tsv'
with src.open(encoding='utf-8-sig') as f:old=list(csv.DictReader(f,delimiter='\t'))
assert len(old)==461
out=[]
for w in json.loads((N/'next_queue.json').read_text()):
 if w['id'] not in ['bg1k_1059','bg1k_1266','bg1k_1900']:continue
 s=int(w['start0']);e=int(w['end0']);same=[r for r in old if r['chrom']==w['chrom']];over=[r for r in same if int(r['start'])<e and int(r['end'])>s]
 distance=lambda r:max(0,int(r['start'])-e,s-int(r['end']))
 nearest=min(same,key=distance) if same else None
 out.append(dict(id=w['id'],chrom=w['chrom'],start0=s,end0=e,original_region_overlaps=[dict(chrom=r['chrom'],start=r['start'],end=r['end'],status=r['integrated_review_status']) for r in over],nearest_original_region=None if nearest is None else dict(start=nearest['start'],end=nearest['end'],distance_bp=distance(nearest)),legacy_wide_ATAC_reusable=bool(over)))
report=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),original_region_count=len(old),results=out,limitations=['No overlap establishes spatial independence from this 461-region set only.','Existing regional ATAC scores cannot be assigned to different loci.','No original veto has been removed or reversed.'])
(N/'legacy_region_overlap_review.json').write_text(json.dumps(report,indent=2));print(json.dumps(out,indent=2))
