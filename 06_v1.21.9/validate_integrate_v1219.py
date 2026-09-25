from pathlib import Path
import csv,json,numpy as np
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.9'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
assert all(r['returncode']==0 for r in json.loads((O/'conservation_execution.json').read_text()))
rows=read(O/'local_windows_conservation_rechecked.tsv');assert len(rows)==2091
key=lambda r:(r['chrom'],r['start'],r['end'])
assert len({key(r) for r in rows})==2091
selected=[r for r in rows if r['selected_1kb_previously']=='True'];assert len(selected)==3
for r in selected:
 a=np.load(O/(r['genes'].split('/')[0]+'.conservation_ROS.npz'));s=int(r['start'])-int(a['start']);v=a['published'][s:s+1000].tolist()
 assert all(np.isfinite(v))
 direct=max(sum(v[i:i+50])/50 for i in range(951))
 assert abs(direct-float(r['published_phyloP_max50bp_mean']))<1e-9
cells=read(P/'06_v1.21.7/cell_marker_masks.tsv');qc={r['barcode'] for r in cells if r['all_QC']=='1'}
assign=read(O/'RNA_cluster_assignments.tsv');assert len(assign)==len(qc)==5849;assert {r['barcode'] for r in assign}==qc
cols=[c for c in assign[0] if c.endswith('_T_enriched')]
groups={c:{r['barcode'] for r in assign if r[c]=='1'} for c in cols};inter=set.intersection(*groups.values());union=set.union(*groups.values());assert(len(inter),len(union))==(155,157)
targets=json.loads((P/'06_v1.21.7/fragment_targets.json').read_text());quant=read(O/'RNA_cluster_local_fragment_counts.tsv');masksets={'T_enriched_all_four_runs':inter,'T_enriched_any_run':union,**{k.removesuffix('_T_enriched'):v for k,v in groups.items()}}
for t in targets:
 fragments=[line.split('\t') for line in (P/'06_v1.21.7'/(t['name'].replace('/','_')+'.fragments.tsv')).read_text().splitlines() if line]
 assert len({tuple(r[:4]) for r in fragments})==len(fragments)
 for q in [q for q in quant if q['target']==t['name']]:
  hits=[r for r in fragments if r[3] in masksets[q['mask']] and t['start']<=int(r[1])<t['end']]
  assert len(hits)==int(q['unique_fragment_starts'])
  assert len({r[3] for r in hits})==int(q['nuclei_with_fragment_start'])
sv=read(O/'NPNT_span_genotypes.tsv');review=json.loads((O/'NPNT_span_review.json').read_text())
for s in review['summary']:
 inside=[r for r in sv if r['sample']==s['sample'] and r['inside']=='True'];hq=[r for r in inside if r['high_quality_het']=='True']
 assert len(inside)==s['inside_sites'];assert len(hq)==s['inside_high_quality_hets'];assert len({r['interval'] for r in hq})==s['inside_intervals_with_high_quality_het']
 assert all(int(r['DP'])>=8 and int(r['GQ'])>=20 and len(set(r['GT'].split('/')))>1 for r in hq)
for r in rows:
 r['legacy_phylop_columns_model']='pilot; use published_phyloP_* for updated calibration'
 if r['selected_1kb_previously']=='True':
  r['T_cell_accessibility_status']='no_fragment_starts_in_155_consensus_or_157_union_RNA_T_enriched_nuclei; not_proof_of_closed_chromatin'
 r['NPNT_regional_SV_context']='catalogue_DEL_challenged_by_distributed_heterozygosity; unresolved_read_level' if r['genes'].startswith('NPNT/') else 'not_applicable'
 r['previous_SV_review']=r['SV_review']
 if r['genes'].startswith('NPNT/') and r['selected_1kb_previously']=='True':r['SV_review']='catalogue_DEL_challenged_by_distributed_heterozygosity; unresolved_read_level'
write(O/'local_windows_integrated_v1219.tsv',rows)
front=[r for r in rows if r['exploratory_Pareto_front']=='True'];assert len(front)==64
write(O/'local_tradeoff_frontier_annotated_v1219.tsv',front)
stats=[]
for g in sorted({r['genes'] for r in front}):
 rr=[r for r in front if r['genes']==g];v=[float(r['published_phyloP_max50bp_mean']) for r in rr]
 stats.append(dict(genes=g,windows=len(rr),published_max50_min=min(v),published_max50_max=max(v)))
result=dict(status='passed',checks=['Three phyloP executions successful; original/new WIG position sets checked during projection','2091 unique local windows; selected 1kb rolling means independently recomputed with Python sums','5849 RNA assignments match QC barcode set exactly; T intersection155/union157','All 42 fragment summaries recomputed from cached rows; no duplicate coordinate/barcode records','NPNT sample counts and high-quality heterozygote criteria recomputed','64 previous Pareto windows preserved; conservation annotated without new threshold or rank'],frontier=stats,limitations=['Checks validate derived arithmetic and bookkeeping, not source experimental quality','Fragment query completeness retains v1.21.7 assumptions','No all-461 rescoring with the published conservation model'])
(O/'integration_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
