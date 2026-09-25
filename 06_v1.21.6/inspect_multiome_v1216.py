from pathlib import Path
import h5py,numpy as np,json,collections,csv
O=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/06_v1.21.6')
with h5py.File(O/'GSE244116_filtered_feature_bc_matrix.h5') as f:
 m=f['matrix'];features=m['features'];types=features['feature_type'].asstr()[:];names=features['name'].asstr()[:];ids=features['id'].asstr()[:];barcodes=m['barcodes'].asstr()[:];shape=m['shape'][:]
 markers=['CD3D','CD3E','CD3G','TRAC','CD247','CD4','CD8A','CD8B','PTPRC','MS4A1','CD79A','LYZ','NKG7']
 index={i:n for i,n in enumerate(names) if n in markers and types[i]=='Gene Expression'}
 counts=np.zeros((len(index),len(barcodes)),dtype=int);lookup={i:j for j,i in enumerate(index)};indptr=m['indptr'][:];indices=m['indices'][:];data=m['data'][:]
 for c in range(len(barcodes)):
  for i,v in zip(indices[indptr[c]:indptr[c+1]],data[indptr[c]:indptr[c+1]]):
   if i in lookup:counts[lookup[i],c]+=v
 marker_names=list(index.values());marker_rows=[dict(marker=n,nuclei_nonzero=int((counts[j]>0).sum()),total_UMI=int(counts[j].sum())) for j,n in enumerate(marker_names)]
 with (O/'multiome_marker_detection.tsv').open('w') as out:
  wr=csv.DictWriter(out,fieldnames=list(marker_rows[0]),delimiter='\t');wr.writeheader();wr.writerows(marker_rows)
 core=[j for j,n in enumerate(marker_names) if n in ['CD3D','CD3E','CD3G','TRAC']];support=(counts[core]>0).sum(axis=0)
 result=dict(shape=shape.tolist(),nuclei=len(barcodes),feature_types=dict(collections.Counter(types)),feature_metadata_keys=list(features),marker_detection=marker_rows,nuclei_detecting_at_least_two_core_T_markers=int((support>=2).sum()),limitation='Marker detection inventory only; no validated cell annotation, no locus accessibility measurement, and no normal T-cell or CAR-T validation. Peak matrix cannot establish absence of accessibility at loci outside called peaks.')
 (O/'multiome_inventory.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
