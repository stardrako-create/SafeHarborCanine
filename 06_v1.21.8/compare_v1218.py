from pathlib import Path
import csv,json,collections
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.21.8'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
old=read(P/'06_v1.21.5/candidates_scored_coherent76.tsv');new=read(O/'candidates_scored_v1218.tsv');key=lambda r:(r['chrom'],r['start'],r['end']);prior={key(r):r for r in old};changes=[]
for r in new:
 before=prior[key(r)]
 for k,v in r.items():
  if v!=before[k]:changes.append(dict(chrom=r['chrom'],start=r['start'],end=r['end'],genes=r['left_gene']+'/'+r['right_gene'],field=k,previous=before[k],current=v))
with (O/'scoring_changes.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(changes[0]),delimiter='\t');wr.writeheader();wr.writerows(changes)
epic=read(P/'06_v1.21.5/epic_nonquiescent_ROS_blocks.tsv');annotations=collections.defaultdict(list)
for r in epic:
 if 1<=int(r['state'])<=7:annotations[r['genes']].append((int(r['ros_start']),int(r['ros_end'])))
def overlap(iv,a,b):
 end=-1;total=0
 for s,e in sorted((max(a,s),min(b,e)) for s,e in iv if s<b and e>a):total+=max(0,e-max(s,end));end=max(end,e)
 return total
for r in new:
 g=r['left_gene']+'/'+r['right_gene'];assessed=g in annotations;r['epic_regional_review']='mapped_blocks_reviewed' if assessed else 'not_assessed';r['epic_regional_promoter_enhancer_bp']=str(overlap(annotations[g],int(r['start']),int(r['end']))) if assessed else ''
 r['integrated_review_status']=('excluded_by_legacy_checks' if r['evaluation_status']=='excluded' else 'regulatory_conflict_requires_local_resolution' if assessed and int(r['epic_regional_promoter_enhancer_bp'])>0 else 'insufficient_evidence')
 r['functional_validation']='not_validated';r['scope_note']='Regional annotation conflict does not exclude every possible subwindow.'
with (O/'candidates_integrated_review_v1218.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(new[0]),delimiter='\t');wr.writeheader();wr.writerows(new)
summary=dict(candidates=len(new),legacy_status_counts=dict(collections.Counter(r['evaluation_status'] for r in new)),status_changes=sum(c['field']=='evaluation_status' for c in changes),hard_veto_changes=sum(c['field']=='hard_veto' for c in changes),changed_fields=dict(collections.Counter(c['field'] for c in changes)),integrated_status_counts=dict(collections.Counter(r['integrated_review_status'] for r in new)),legacy_survivors=[{k:r[k] for k in ['left_gene','right_gene','final_score','atac_mean_percentile','atac_mean_narrow_percentile','epic_regional_promoter_enhancer_bp','integrated_review_status']} for r in new if r['evaluation_status']=='passes_recorded_checks'],validated_safe_harbors=0)
(O/'scoring_comparison_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
