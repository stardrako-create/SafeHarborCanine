from pathlib import Path
import csv,json,h5py,numpy as np,networkx as nx,sklearn
from scipy.sparse import csc_matrix,diags
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.9'
cells={r['barcode']:r for r in csv.DictReader((P/'06_v1.21.7/cell_marker_masks.tsv').open(),delimiter='\t')}
with h5py.File(P/'06_v1.21.6/GSE244116_filtered_feature_bc_matrix.h5') as f:
 m=f['matrix'];names=m['features/name'].asstr()[:];types=m['features/feature_type'].asstr()[:];bc=m['barcodes'].asstr()[:];matrix=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
keep=np.array([cells[b]['all_QC']=='1' for b in bc]);bc=bc[keep];gn=names[types=='Gene Expression'];raw=matrix[types=='Gene Expression'][:,keep].tocsr();del matrix
counts=np.asarray(raw.sum(0)).ravel();norm=(raw@diags(10000/counts)).tocsr();norm.data=np.log1p(norm.data);mean=np.asarray(norm.mean(1)).ravel();var=np.asarray(norm.multiply(norm).mean(1)).ravel()-mean**2;det=np.asarray(raw.getnnz(axis=1)).ravel()
# Normalize variance within mean-expression bins; annotate method explicitly rather than calling this Seurat VST.
eligible=(det>=10)&np.array([not n.startswith(('MT-','mt-','RPS','RPL')) for n in gn]);score=np.full(len(gn),-np.inf);ordered=np.flatnonzero(eligible)[np.argsort(mean[eligible])]
for ix in np.array_split(ordered,20):
 vv=np.log(np.maximum(var[ix],1e-12));score[ix]=(vv-vv.mean())/(vv.std()+1e-12)
hvg=np.argsort(score)[-2000:];x=norm[hvg].toarray().T.astype(np.float32);x=(x-x.mean(0))/(x.std(0)+1e-6);np.clip(x,-10,10,out=x)
with threadpool_limits(limits=2):
 pca=PCA(n_components=30,svd_solver='randomized',random_state=42);pcs=pca.fit_transform(x);dist,neighbors=NearestNeighbors(n_neighbors=21).fit(pcs).kneighbors(pcs)
graph=nx.Graph();graph.add_nodes_from(range(len(bc)));scale=float(np.median(dist[:,1:]));assert scale>0
for i in range(len(bc)):
 for d,j in zip(dist[i,1:],neighbors[i,1:]):graph.add_edge(i,int(j),weight=float(np.exp(-d/scale)))
def expr(n):return np.asarray(raw[gn==n].sum(0)).ravel()
core=np.sum([expr(n)>0 for n in ['CD3D','CD3E','CD3G']],axis=0);support=np.sum([expr(n)>0 for n in ['CD247','LCK','CD2']],axis=0);my=np.sum([expr(n)>0 for n in ['LYZ','LST1','CSF1R']],axis=0);bcell=np.sum([expr(n)>0 for n in ['MS4A1','CD79A','CD79B']],axis=0)
labels={};selected={};profiles=[];marker_profiles=[]
for resolution in [.5,1.]:
 for seed in [0,42]:
  label=f'res{resolution}_seed{seed}';communities=nx.community.louvain_communities(graph,resolution=resolution,seed=seed);communities=sorted(communities,key=lambda c:min(c));lab=np.zeros(len(bc),dtype=int);Tset=set()
  for k,c in enumerate(communities):
   ix=np.array(sorted(c));lab[ix]=k;cf=float(np.mean(core[ix]>0));sf=float(np.mean(support[ix]>0));mf=float(np.mean(my[ix]>=2));bf=float(np.mean(bcell[ix]>=2));Tlike=len(ix)>=20 and cf>=.3 and sf>=.3 and mf<.2 and bf<.2
   if Tlike:Tset.update(ix.tolist())
   profiles.append(dict(run=label,cluster=k,nuclei=len(ix),core_T_fraction=cf,support_T_fraction=sf,myeloid_2plus_fraction=mf,B_2plus_fraction=bf,T_enriched_by_rule=Tlike))
   for n in ['CD3D','CD3E','CD3G','CD247','LCK','CD2','CD4','CD8A','CD8B','LYZ','LST1','CSF1R','MS4A1','CD79A','NKG7','MKI67']:
    v=expr(n)[ix];marker_profiles.append(dict(run=label,cluster=k,gene=n,detected_fraction=float(np.mean(v>0)),mean_UMI=float(v.mean())))
  labels[label]=lab;selected[label]=Tset;print(label,'clusters',len(communities),'T-enriched',len(Tset),flush=True)
def write(path,rows):
 with path.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write(O/'RNA_cluster_profiles.tsv',profiles);write(O/'RNA_cluster_marker_profiles.tsv',marker_profiles)
write(O/'RNA_cluster_assignments.tsv',[dict(barcode=b,**{k:int(v[i]) for k,v in labels.items()},**{k+'_T_enriched':int(i in s) for k,s in selected.items()}) for i,b in enumerate(bc)])
common=set.intersection(*selected.values());union=set.union(*selected.values());quant=[]
for mask,ix in [('T_enriched_all_four_runs',common),('T_enriched_any_run',union)]+list(selected.items()):
 barcodes={bc[i] for i in ix}
 for target in json.loads((P/'06_v1.21.7/fragment_targets.json').read_text()):
  path=P/'06_v1.21.7'/(target['name'].replace('/','_')+'.fragments.tsv');rs=[l.split('\t') for l in path.read_text().splitlines() if l];starts=[r for r in rs if r[3] in barcodes and target['start']<=int(r[1])<target['end']]
  quant.append(dict(mask=mask,nuclei=len(barcodes),target=target['name'],unique_fragment_starts=len(starts),nuclei_with_fragment_start=len({r[3] for r in starts})))
write(O/'RNA_cluster_local_fragment_counts.tsv',quant)
summary=dict(nuclei=len(bc),genes=len(gn),hvg_count=len(hvg),pca_components=30,knn=20,versions=dict(networkx=nx.__version__,sklearn=sklearn.__version__),hvg_method='log1p library-size-normalized RNA; variance z-score within 20 mean-expression bins; not Seurat VST',clustering='RNA-only weighted kNN Louvain; not authors multimodal WNN reproduction',T_enriched_rule='>=20 nuclei; >=30% expressing >=1 CD3D/E/G and >=30% expressing >=1 CD247/LCK/CD2; <20% with >=2 myeloid or B markers. Exploratory cluster label, not validated annotation.',intersection_n=len(common),union_n=len(union),run_T_counts={k:len(v) for k,v in selected.items()},ARI={a+' vs '+b:adjusted_rand_score(labels[a],labels[b]) for i,a in enumerate(labels) for b in list(labels)[i+1:]},limitations=['Single tumor donor','No independent barcode-level ground truth','No doublet detection or full ATAC QC','Cluster annotation thresholds exploratory; not purity estimates','Target ATAC counts not used for clustering or annotation'])
(O/'RNA_cluster_summary.json').write_text(json.dumps(summary,indent=2));np.savez_compressed(O/'RNA_cluster_embedding.npz',barcodes=bc.astype(str),pcs=pcs,hvg_names=gn[hvg].astype(str));print(json.dumps(summary,indent=2))
