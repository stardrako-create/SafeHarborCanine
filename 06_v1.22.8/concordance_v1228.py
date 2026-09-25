from pathlib import Path
import csv,gzip,json,collections,hashlib,itertools
import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import linear_sum_assignment
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.8';O.mkdir(exist_ok=True)
with gzip.open(P/'06_v1.22.3/GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz','rt') as f:published=list(csv.DictReader(f,delimiter='\t'))
with (P/'06_v1.22.7/gene_counts_unstranded.tsv').open() as f:r=csv.DictReader(f,delimiter='\t');names=r.fieldnames[1:];raw=list(r)
assert set(names)==set(published[0])-{'Gene'}
symbols=collections.defaultdict(list)
for l in (P/'06_v1.22.7/STAR_index/geneInfo.tab').read_text().splitlines()[1:]:
 v=l.split('\t');symbols[v[1]].append(v[0])
pubsymbols=collections.Counter(r['Gene'] for r in published);rawlookup={r['gene_id']:r for r in raw};pairs=[]
for r in published:
 g=r['Gene']
 if pubsymbols[g]==1 and len(symbols[g])==1 and symbols[g][0] in rawlookup:pairs.append((g,symbols[g][0],r,rawlookup[symbols[g][0]]))
X=np.array([[float(r[n]) for n in names] for g,i,p,r in pairs]);Y=np.array([[float(p[n]) for n in names] for g,i,p,r in pairs]);totals=np.array([sum(float(r[n]) for r in raw) for n in names]);X=X/totals*1e6
assert np.isfinite(X).all() and np.isfinite(Y).all() and (X>=0).all() and (Y>=0).all()
def write(name,rows):
 with (O/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write('exact_symbol_mapping.tsv',[dict(symbol=g,ROS_gene_id=i) for g,i,p,r in pairs])
results=[];assess=[]
for label,mask in [('all_exact',np.ones(len(pairs),dtype=bool)),('CPM1_both_at_least3_samples',((X>=1).sum(axis=1)>=3)&((Y>=1).sum(axis=1)>=3)),('CPM5_both_at_least3_samples',((X>=5).sum(axis=1)>=3)&((Y>=5).sum(axis=1)>=3))]:
 x=X[mask];y=Y[mask];lx=np.log2(x+1);ly=np.log2(y+1);cx=lx-lx.mean(axis=1,keepdims=True);cy=ly-ly.mean(axis=1,keepdims=True)
 for metric in ['spearman_abundance','pearson_gene_centered_log2_CPM_plus1']:
  matrix=np.array([[spearmanr(x[:,i],y[:,j]).statistic if metric=='spearman_abundance' else np.corrcoef(cx[:,i],cy[:,j])[0,1] for j in range(9)] for i in range(9)])
  assert np.isfinite(matrix).all()
  for i in range(9):
   for j in range(9):results.append(dict(gene_set=label,genes=int(mask.sum()),metric=metric,reprocessed_sample=names[i],published_sample=names[j],correlation=float(matrix[i,j])))
  ri,ci=linear_sum_assignment(-matrix)
  assess.append(dict(gene_set=label,genes=int(mask.sum()),metric=metric,diagonal_min=float(np.diag(matrix).min()),diagonal_max=float(np.diag(matrix).max()),diagonal_best_count=int(sum(matrix[i,i]>=max(matrix[i]) for i in range(9))),best_assignment=[dict(reprocessed=names[i],published=names[j]) for i,j in zip(ri,ci)],identity_assignment=bool(np.all(ri==ci)),diagonal_sum=float(np.trace(matrix)),optimal_sum=float(matrix[ri,ci].sum())))
write('cross_sample_concordance.tsv',results)
unmapped=[dict(published_symbol=r['Gene'],reason='missing_exact_ROS_symbol' if not symbols[r['Gene']] else 'ambiguous_symbol') for r in published if not (pubsymbols[r['Gene']]==1 and len(symbols[r['Gene']])==1 and symbols[r['Gene']][0] in rawlookup)]
if unmapped:write('unmapped_published_symbols.tsv',unmapped)
meta=list(csv.DictReader((P/'06_v1.22.4/sample_metadata_audit.tsv').open(),delimiter='\t'));write('metadata_conflicts_preserved.tsv',[r for r in meta if r['metadata_conflict']=='specific_source_or_tissue_conflicts'])
summary=dict(exact_unique_genes=len(pairs),published_genes=len(published),unmapped=len(unmapped),assessments=assess,limitations=['Exact unique gene symbols only; not proof of one-to-one orthologous exon definitions','Same sequencing dataset: concordance is provenance consistency, not independent biological validation','Cannot resolve physical donor identity or mislabeled samples shared by raw and published files','log2(CPM+1) only for descriptive similarity; no inferential ratios or differential expression','Sensitivity gene filters are descriptive, not candidate filters','No labels reassigned; four metadata conflicts preserved'])
(O/'concordance_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
