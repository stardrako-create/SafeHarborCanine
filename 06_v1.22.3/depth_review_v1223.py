from pathlib import Path
import csv,json,numpy as np,h5py
from scipy.sparse import csc_matrix
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.3';O.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
cells={r['barcode']:r for r in read(P/'06_v1.21.7/cell_marker_masks.tsv')};assign=read(P/'06_v1.21.9/RNA_cluster_assignments.tsv');cols=[c for c in assign[0] if c.endswith('_T_enriched')];groups=[{r['barcode'] for r in assign if r[c]=='1'} for c in cols];cons=set.intersection(*groups);union=set.union(*groups);qc={r['barcode'] for r in assign};other=qc-union
metrics=['RNA_features','RNA_UMI','ATAC_peak_features','ATAC_peak_counts'];rows=[]
for name,bars in [('T_consensus',cons),('T_union',union),('other_QC_excluding_T_union',other),('all_QC',qc)]:
 for metric in metrics:
  v=np.array([int(cells[b][metric]) for b in bars]);rows.append(dict(group=name,nuclei=len(bars),metric=metric,total=int(v.sum()),minimum=int(v.min()),p10=float(np.percentile(v,10)),q25=float(np.percentile(v,25)),median=float(np.median(v)),q75=float(np.percentile(v,75)),p90=float(np.percentile(v,90)),maximum=int(v.max())))
# Recompute published count/feature metrics from H5, ensuring no cached-column drift.
with h5py.File(P/'06_v1.21.6/GSE244116_filtered_feature_bc_matrix.h5') as f:
 m=f['matrix'];types=m['features/feature_type'].asstr()[:];bc=m['barcodes'].asstr()[:];mat=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
 for typ,prefix,countname in [('Gene Expression','RNA','RNA_UMI'),('Peaks','ATAC_peak','ATAC_peak_counts')]:
  data=mat[types==typ];counts=np.asarray(data.sum(0)).ravel();features=np.asarray(data.getnnz(0)).ravel()
  assert all(int(counts[i])==int(cells[b][countname]) and int(features[i])==int(cells[b][prefix+'_features']) for i,b in enumerate(bc))
# Compare local starts with cells of similar peak-matrix depth. Conditional resampling only,
# not a biological significance test and not full genomic fragment normalization.
barcodes=sorted(qc);index={b:i for i,b in enumerate(barcodes)};depth=np.array([int(cells[b]['ATAC_peak_counts']) for b in barcodes]);Tidx=np.array([index[b] for b in sorted(cons)]);Cidx=np.array([index[b] for b in sorted(other)])
edges=np.quantile(depth,np.linspace(0,1,11));bins=np.searchsorted(edges[1:-1],depth,side='right');targets=[]
for r in read(P/'06_v1.22.0/frontier_T_evidence_v1220.tsv'):
 if r['canfam6_mapping_status']=='two_routes_unique_identical_contiguous':targets.append(dict(name=r['genes']+':'+r['start']+'-'+r['end'],genes=r['genes'],kind='frontier_full_1kb',start=int(r['canfam6_start']),end=int(r['canfam6_end']),file=P/'06_v1.22.0'/(r['genes'].split('/')[0]+'_frontier_span.fragments.tsv')))
for r in json.loads((P/'06_v1.21.7/fragment_targets.json').read_text()):targets.append(dict(name=r['name'],genes=r['name'],kind=r['kind'],start=r['start'],end=r['end'],file=P/'06_v1.21.7'/(r['name'].replace('/','_')+'.fragments.tsv')))
counts=np.zeros((len(barcodes),len(targets)),dtype=np.int32)
for j,t in enumerate(targets):
 for line in t['file'].read_text().splitlines():
  if not line:continue
  r=line.split('\t')
  if r[3] in index and t['start']<=int(r[1])<t['end']:counts[index[r[3]],j]+=1
rng=np.random.default_rng(42);R=2000;samplecounts=np.zeros((R,len(targets)),dtype=np.int32);sampledepth=np.zeros(R);strata=[]
for k in range(10):
 need=sum(bins[Tidx]==k);pool=Cidx[bins[Cidx]==k];assert len(pool)>=need;strata.append(dict(bin=k,T_nuclei=int(need),control_pool=int(len(pool))))
for i in range(R):
 selected=np.concatenate([rng.choice(Cidx[bins[Cidx]==k],size=sum(bins[Tidx]==k),replace=False) for k in range(10)]);assert len(selected)==155 and len(set(selected))==155
 samplecounts[i]=counts[selected].sum(0);sampledepth[i]=depth[selected].sum()
local=[]
for j,t in enumerate(targets):
 v=samplecounts[:,j];local.append(dict(target=t['name'],kind=t['kind'],T_nuclei=155,T_fragment_starts=int(counts[Tidx,j].sum()),T_nuclei_with_start=int((counts[Tidx,j]>0).sum()),other_QC_fragment_starts=int(counts[Cidx,j].sum()),matched_other_sample_median=float(np.median(v)),matched_other_sample_p025=float(np.percentile(v,2.5)),matched_other_sample_p975=float(np.percentile(v,97.5)),fraction_matched_samples_with_zero=float(np.mean(v==0))))
for name,rr in [('depth_group_summary.tsv',rows),('depth_matched_local_counts.tsv',local)]:
 with (O/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
summary=dict(QC_nuclei=len(qc),T_consensus=155,T_union=157,other_QC=len(other),T_fraction_nuclei=len(cons)/len(qc),T_fraction_peak_matrix_counts=float(depth[Tidx].sum()/depth.sum()),T_peak_matrix_counts=int(depth[Tidx].sum()),matched_sample_peak_count_median=float(np.median(sampledepth)),matched_sample_peak_count_p025=float(np.percentile(sampledepth,2.5)),matched_sample_peak_count_p975=float(np.percentile(sampledepth,97.5)),sampling=dict(seed=42,replicates=R,strata=strata,decile_boundaries=edges.tolist()),cached_metrics_recomputed_from_H5='passed',local_targets=len(targets),limitations=['Peak-matrix count is a depth proxy, not total genomic fragments or FRiP','Count/feature QC does not assess TSS enrichment, nucleosome signal, doublets or barcode identity','Depth-decile conditional resampling is descriptive; intervals are distributions of resampled other cells, not confidence intervals or biological significance','Controls exclude union T cells but include heterogeneous tumor/microenvironment populations','One tumor donor; no inference of normal T accessibility or CAR-T safety','Only37 full contiguous frontier windows use this comparison;27 partially mapped alternatives not directly comparable1kb'])
(O/'depth_review_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(json.dumps([r for r in rows if r['metric']=='ATAC_peak_counts'],indent=2));print(json.dumps([r for r in local if r['kind']!='frontier_full_1kb'],indent=2))
