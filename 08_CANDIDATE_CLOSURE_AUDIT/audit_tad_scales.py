from pathlib import Path
import csv,json,hashlib,shutil
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT'
src=P/'04_tracks_processadas/ROS_Cfam_1.0/zenodo_upload/HiC_tad_insulation.tsv'
rows=list(csv.DictReader(src.open(),delimiter='\t'))
cols=[k for k in rows[0] if k.startswith('is_boundary_')]
risk={r['gene_symbol'] for r in csv.DictReader((P/'05_SHIP/canine_risk_genes.tsv').open(),delimiter='\t')}
genes=[l.rstrip().split('\t') for l in (P/'05_SHIP/canine_all_genes.bed').open()]
windows=[r for r in csv.DictReader((P/'07_FINAL_CANDIDATES_2026-09-19/shortlist_evidence.tsv').open(),delimiter='\t') if r['window_id'] in ('w01','w11')]
out=[]
for w in windows:
 c,s,e=w['chrom'],int(w['start']),int(w['end'])
 for scale in ['union']+cols:
  bounds=sorted((int(r['start']),int(r['end'])) for r in rows if r['chrom']==c and any(r[k].lower()=='true' for k in (cols if scale=='union' else [scale])))
  merged=[]
  for a,b in bounds:
   if merged and a<=merged[-1][1]:merged[-1][1]=max(b,merged[-1][1])
   else:merged.append([a,b])
  intervals=[(a[1],b[0]) for a,b in zip(merged,merged[1:]) if a[1]<b[0]]
  own=[(a,b) for a,b in intervals if a<=s and e<=b]
  assert len(own)<=1
  hits=sorted({g[3] for g in genes if own and g[0]==c and int(g[1])<own[0][1] and int(g[2])>own[0][0] and g[3] in risk})
  out.append(dict(window=w['window_id'],scale=scale,chrom=c,proxy_start=own[0][0] if own else None,proxy_end=own[0][1] if own else None,risk_genes=hits,status='risk_catalogue_overlap' if hits else 'no_catalogue_hit_in_boundary_gap' if own else 'no_full_window_assignment'))
result={'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'results':out,'limitations':['Boundary gaps are proxies, not validated domains or contact maps in canine T cells.','Union of scales can split a larger-scale boundary gap and hide risk genes detected at another scale.','Risk catalogue overlap does not establish a regulatory interaction; no overlap does not establish safety.']}
(O/'tad_scale_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(__file__,O/Path(__file__).name)
print(json.dumps(result,indent=2))
