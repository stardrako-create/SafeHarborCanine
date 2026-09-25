from pathlib import Path
import csv,json,hashlib,collections,datetime
P=Path(__file__).parent.parent;O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922'
src=P/'04_tracks_processadas/ROS_Cfam_1.0/zenodo_upload/HiC_tad_insulation.tsv';rows=list(csv.DictReader(src.open(),delimiter='\t'));cols=[k for k in rows[0] if k.startswith('is_boundary_')]
risk={r['gene_symbol'] for r in csv.DictReader((P/'05_SHIP/canine_risk_genes.tsv').open(),delimiter='\t')};genes=[l.rstrip().split('\t') for l in (P/'05_SHIP/canine_all_genes.bed').open()]
windows=list(csv.DictReader((O/'transcript_and_boundary_review.tsv').open(),delimiter='\t'));out=[]
for w in windows:
 c,s,e=w['chrom'],int(w['start0']),int(w['end0'])
 for scale in ['union']+cols:
  bounds=sorted((int(r['start']),int(r['end'])) for r in rows if r['chrom']==c and any(r[k].lower()=='true' for k in (cols if scale=='union' else [scale])))
  merged=[]
  for a,b in bounds:
   if merged and a<=merged[-1][1]:merged[-1][1]=max(b,merged[-1][1])
   else:merged.append([a,b])
  own=[(a[1],b[0]) for a,b in zip(merged,merged[1:]) if a[1]<=s and e<=b[0]];assert len(own)<=1
  hits=sorted({g[3] for g in genes if own and g[0]==c and int(g[1])<own[0][1] and int(g[2])>own[0][0] and g[3] in risk})
  out.append(dict(id=w['id'],group=w['group'],scale=scale,proxy_start=own[0][0] if own else None,proxy_end=own[0][1] if own else None,risk_genes=hits,status='risk_catalogue_overlap' if hits else 'no_catalogue_hit_in_boundary_gap' if own else 'no_full_window_assignment'))
by={w['id']:[r for r in out if r['id']==w['id']] for w in windows};hitids=sorted(k for k,v in by.items() if any(r['risk_genes'] for r in v));missing=sorted(k for k,v in by.items() if any(r['status']=='no_full_window_assignment' for r in v));clear=sorted(k for k,v in by.items() if all(r['status']=='no_catalogue_hit_in_boundary_gap' for r in v))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),reviewed=len(windows),risk_hit_ids=hitids,missing_assignment_ids=missing,no_hit_all_scales_ids=clear,results=out,limitations=['Same boundary-gap proxy method as existing w01/w11 audit; not T-cell contact validation','Risk and missing-assignment categories can overlap','Union can hide risk at coarser scales; use all scale results','No hit is not proof of regulatory isolation'])
(O/'TAD_multiscale_review.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='results'},indent=2))
