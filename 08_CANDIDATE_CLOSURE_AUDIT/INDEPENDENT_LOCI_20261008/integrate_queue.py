from pathlib import Path
import json,hashlib,datetime
N=Path(__file__).resolve().parent
files=['TAD_clear_queue.json','EpiC_roundtrip_and_state_review.json','UU_projection/projection_review.json','SV_all_windows/status.json','sequence_RRBS_review.json','variant_catalogue_query_drop_samples/status.json']
d={f:json.loads((N/f).read_text(encoding='utf-8-sig')) for f in files}
e={r['id']:r for r in d[files[1]]['state_review']};u={r['id']:r for r in d[files[2]]['results']};rr={r['id']:r for r in d[files[4]]['results']}
rows=[]
for r in d[files[0]]:
 k=f"{r['chrom']}_{r['start0']}_{r['end0']}";h=[v for v in d[files[3]]['hits'] if v['query']==k];ep=e.get(k)
 clear=ep is not None and ep['roundtrip_pass'] and not ep['promoter_enhancer_tissues']
 rows.append(dict(**r,coordinate_id=k,SV_hits=h,EpiC_clear_with_roundtrip=clear,UU_exact=u[k]['status']=='eligible_exact_window_for_variant_query',RRBS_observed_bp=rr[r['id']]['RRBS_observed_bp'],RRBS_observed_mean=rr[r['id']]['RRBS_observed_mean'],softmasked_bp=rr[r['id']]['softmasked_bp'],next_priority=clear and not h,approved=False))
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':{f:hashlib.sha256((N/f).read_bytes()).hexdigest() for f in files},'counts':{'windows':len(rows),'with_SV_hit':sum(bool(r['SV_hits']) for r in rows),'EpiC_clear_no_SV_hit':sum(r['next_priority'] for r in rows),'exact_UU_EpiC_clear_no_SV_hit':sum(r['next_priority'] and r['UU_exact'] for r in rows),'completed_variant_queries':len(d[files[5]]['completed']),'query_errors':len(d[files[5]]['errors'])},'results':rows,'limitations':['No catalogue hit does not establish absence of variation.','ATAC narrow and wide legacy gates, conservation and specificity remain unverified for new loci.','Sparse or absent RRBS cannot demonstrate low methylation.','This is prioritization, not computational approval.']}
(N/'integrated_queue_20261008.json').write_text(json.dumps(out,indent=2))
(N/'next_queue.json').write_text(json.dumps([r for r in rows if r['next_priority']],indent=2))
print(json.dumps(out['counts']))
print(json.dumps([{k:r[k] for k in ['id','coordinate_id','atac_percentile','UU_exact','RRBS_observed_bp','RRBS_observed_mean']} for r in rows if r['next_priority']],indent=2))
