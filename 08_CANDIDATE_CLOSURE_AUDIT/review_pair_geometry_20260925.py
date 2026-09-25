from pathlib import Path
import json,collections,datetime,hashlib
A=Path(__file__).resolve().parent/'SV_PBGV000010/alignment_chunked';I=A/'joint_alternative_panel_20260925';O=A/'pair_geometry_20260925';O.mkdir(exist_ok=False)
hits=json.loads((I/'all_hits.json').read_text());meta=json.loads((I/'genomic_coordinates.json').read_text());single=json.loads((I/'review.json').read_text());selected={r['molecule'] for r in single['read_results'] if r['best_spanning_AS'] is not None}
# Enumerate reported alignments, collapsing identical genomic placements from overlapping context.
def placements(name,kind):
 out={}
 for h in hits[name]:
  if h['genomic']!=(kind=='genomic'):continue
  if kind=='genomic':m=meta[h['target']];target=m['chrom'];s=m['start0']+h['start0'];e=m['start0']+h['end0']
  else:target=h['target'];s=h['start0'];e=h['end0']
  k=(target,s,e,h['reverse'],h['cigar']);r=dict(h,coordinate_target=target,coordinate_start0=s,coordinate_end0=e)
  if k not in out or h['AS']>out[k]['AS']:out[k]=r
 return list(out.values())
def geometry(a,b):
 if a['coordinate_target']!=b['coordinate_target']:return None
 left,right=sorted([a,b],key=lambda x:x['coordinate_start0']);inward=not left['reverse'] and right['reverse']
 return {'inward':inward,'outer_span_bp':max(a['coordinate_end0'],b['coordinate_end0'])-min(a['coordinate_start0'],b['coordinate_start0']),'score_sum':a['AS']+b['AS']}
results=[]
for molecule in sorted({n.rsplit('_mate',1)[0] for n in hits}):
 names=[molecule+'_mate1',molecule+'_mate2'];assert all(n in hits for n in names)
 candidates={};counts={}
 for kind in ['genomic','junction']:
  pairs=[]
  for a in placements(names[0],kind):
   for b in placements(names[1],kind):
    g=geometry(a,b)
    if g is None or not g['inward']:continue
    if kind=='junction' and not(a['junction_spanning'] or b['junction_spanning']):continue
    pairs.append(dict(g,mate1=a,mate2=b))
  candidates[kind]=pairs;counts[kind]=len(pairs)
 sensitivity=[]
 for maximum in [500,1000,2000]:
  best={}
  for kind,pairs in candidates.items():
   eligible=[r for r in pairs if r['outer_span_bp']<=maximum];score=max((r['score_sum'] for r in eligible),default=None);top=[r for r in eligible if r['score_sum']==score]
   best[kind]={'best_score_sum':score,'top_pair_count':len(top),'best_pairs':top}
  g=best['genomic']['best_score_sum'];j=best['junction']['best_score_sum'];sensitivity.append({'max_outer_span_bp':maximum,'junction_minus_genomic_score':None if j is None or g is None else j-g,**best})
 results.append({'molecule_name':molecule,'previous_spanning_read_selected':molecule in selected,'inward_pair_counts_all_reported_distances':counts,'sensitivity':sensitivity})
report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'molecule_names':len(results),'previous_spanning_molecule_names':len(selected),'results':results,'input_sha256':{f:hashlib.sha256((I/f).read_bytes()).hexdigest() for f in ['all_hits.json','genomic_coordinates.json','reads.fa']},'limitations':['Same 30 selected reads as prior analysis; no additional independent support','Single-read BWA candidates enumerated as pairs, not paired realignment or calibrated fragment likelihood','500/1000/2000bp are explicit descriptive sensitivity bounds, not new acceptance criteria or estimated library distribution','Panel consists of reported alternatives, not exhaustive whole-genome competitive mapping','AS sums are alignment scores, not posterior probabilities or genotypes; mates grouped by name do not prove molecule independence']}
(O/'review.json').write_text(json.dumps(report,indent=2));print(json.dumps([{'name':r['molecule_name'],'sensitivity':[{'span':s['max_outer_span_bp'],'margin':s['junction_minus_genomic_score'],'junction_spans':sorted(set(p['outer_span_bp'] for p in s['junction']['best_pairs']))} for s in r['sensitivity']]} for r in results if r['previous_spanning_read_selected']],indent=2))
