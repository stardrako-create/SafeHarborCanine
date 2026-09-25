from pathlib import Path
import h5py,numpy as np,csv,json
from scipy.sparse import csc_matrix
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.7'
with h5py.File(P/'06_v1.21.6/GSE244116_filtered_feature_bc_matrix.h5') as f:
 m=f['matrix'];ft=m['features'];names=ft['name'].asstr()[:];types=ft['feature_type'].asstr()[:];intervals=ft['interval'].asstr()[:];bc=m['barcodes'].asstr()[:];mat=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
rna=mat[types=='Gene Expression'];atac=mat[types=='Peaks'];nr=np.asarray(rna.getnnz(axis=0)).ravel();na=np.asarray(atac.getnnz(axis=0)).ravel();cr=np.asarray(rna.sum(axis=0)).ravel();ca=np.asarray(atac.sum(axis=0)).ravel();qc=(nr>100)&(nr<30000)&(cr>50)&(cr<50000)&(na>100)&(na<30000)&(ca>50)&(ca<50000)
def expr(n):
 ix=np.flatnonzero((names==n)&(types=='Gene Expression'));return np.asarray(mat[ix].sum(axis=0)).ravel()
core=sum(expr(n)>0 for n in ['CD3D','CD3E','CD3G']);support=sum(expr(n)>0 for n in ['CD247','CD2','LCK']);myeloid=sum(expr(n)>0 for n in ['LYZ','LST1','CSF1R']);bcell=sum(expr(n)>0 for n in ['MS4A1','CD79A','CD79B'])
masks={'all_QC':qc,'T_marker_2plus':qc&(core>=2),'T_marker_strict':qc&(core>=2)&(support>=1)&(myeloid<2)&(bcell<2),'T_marker_broad':qc&(core>=1)&(support>=1)&(myeloid<2)&(bcell<2)}
rows=[]
for i,b in enumerate(bc):rows.append(dict(barcode=b,RNA_features=int(nr[i]),RNA_UMI=int(cr[i]),ATAC_peak_features=int(na[i]),ATAC_peak_counts=int(ca[i]),CD3_genes_detected=int(core[i]),T_support_genes_detected=int(support[i]),myeloid_genes_detected=int(myeloid[i]),B_genes_detected=int(bcell[i]),**{n:int(v[i]) for n,v in masks.items()}))
with (O/'cell_marker_masks.tsv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
aliases={l.split()[0]:l.split()[3] for l in (O/'canFam6.chromAlias.txt').read_text().splitlines() if l and not l.startswith('#')}
targets=[];peakrows=[]
for r in csv.DictReader((O/'local_canfam6_windows.tsv').open(),delimiter='\t'):
 chrom=aliases[r['chrom']];a=int(r['start']);b=int(r['end']);targets.append(dict(name=r['genes'],chrom=chrom,start=a,end=b,kind='candidate'))
 ix=[]
 for j,iv in enumerate(intervals):
  if types[j]!='Peaks' or not iv.startswith(chrom+':'):continue
  s,e=map(int,iv.split(':')[1].split('-'))
  if s<b and e>a:ix.append(j)
 for n,mask in masks.items():peakrows.append(dict(genes=r['genes'],mask=n,nuclei=int(mask.sum()),overlapping_matrix_peaks=len(ix),counts_in_overlapping_whole_peaks=int(mat[ix][:,mask].sum()),caveat='Whole-peak count, not counts limited to candidate interval; zero overlapping peaks is unassessed local signal'))
# Controls use the feature-associated interval midpoint; do not assert that these are verified TSSs.
for n in ['CD3D','CD3E','CD247','LYZ']:
 j=np.flatnonzero((names==n)&(types=='Gene Expression'));assert len(j)==1
 iv=intervals[j[0]];chrom,span=iv.split(':');a,b=map(int,span.split('-'));mid=(a+b)//2
 targets.append(dict(name=n+'_feature_control',chrom=chrom,start=mid-500,end=mid+500,kind='gene_feature_control_not_verified_TSS'))
(O/'fragment_targets.json').write_text(json.dumps(targets,indent=2))
with (O/'candidate_peak_matrix_counts.tsv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(peakrows[0]),delimiter='\t');wr.writeheader();wr.writerows(peakrows)
(O/'cell_mask_summary.json').write_text(json.dumps(dict(total_nuclei=len(bc),masks={n:int(v.sum()) for n,v in masks.items()},method='Author RNA/ATAC count-feature QC. Exploratory transcript-only marker masks fixed without looking at target ATAC counts. Not reproduced WNN clusters or validated cell labels.',core='CD3D/E/G',support='CD247/CD2/LCK',strict='>=2 core and >=1 support; <2 myeloid and <2 B markers',limitations=['One tumor donor','Marker dropout and doublets remain','QC count thresholds do not replace TSS enrichment, nucleosome signal or doublet assessment','No clinical or functional validation']),indent=2));print({n:int(v.sum()) for n,v in masks.items()});print(peakrows)
