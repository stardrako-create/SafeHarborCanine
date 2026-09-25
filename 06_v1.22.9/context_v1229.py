from pathlib import Path
import csv,json,collections,statistics,math,hashlib
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.9';O.mkdir(exist_ok=True)
previous=json.loads((P/'06_v1.22.4/RNA_context_validation.json').read_text());genes=sorted({r['gene'] for r in csv.DictReader((P/'06_v1.22.4/CAR_T_gene_expression_context.tsv').open(),delimiter='\t')});assert len(genes)==40
symbols=collections.defaultdict(list)
for l in (P/'06_v1.22.7/STAR_index/geneInfo.tab').read_text().splitlines()[1:]:
 v=l.split('\t');symbols[v[1]].append(v[0])
with (P/'06_v1.22.7/gene_counts_unstranded.tsv').open() as f:reader=csv.DictReader(f,delimiter='\t');names=reader.fieldnames[1:];matrix={r['gene_id']:r for r in reader}
totals={n:sum(int(r[n]) for r in matrix.values()) for n in names};mapping=[];long=[];conditions=[];pairs=[];overview=[]
for gene in genes:
 ids=symbols[gene];status='exact_unique' if len(ids)==1 and ids[0] in matrix else 'missing_annotation_symbol' if not ids else 'ambiguous_symbol';mapping.append(dict(gene=gene,ROS_gene_ids=';'.join(ids),mapping_status=status,absent_published_CPM=gene in previous['genes_absent']))
 counts={n:int(matrix[ids[0]][n]) for n in names} if status=='exact_unique' else {};cpm={n:counts[n]*1e6/totals[n] for n in counts}
 for n in names:long.append(dict(gene=gene,sample=n,reported_donor=n.split('_')[0],mapping_status=status,raw_count=counts.get(n,''),CPM=cpm.get(n,''),observation='annotation_unresolved' if n not in counts else 'zero_assigned_pairs' if counts[n]==0 else 'positive_assigned_pairs'))
 for label,suffix in [('control','T_cell'),('single_CAR','B7H3_CAR_T'),('dual_CAR','BC_CAR_T')]:
  vs=[cpm[d+'_'+suffix] for d in ['B','E','M']] if counts else [];conditions.append(dict(gene=gene,condition=label,n_reported_donors=len(vs),median_CPM=statistics.median(vs) if vs else '',min_CPM=min(vs) if vs else '',max_CPM=max(vs) if vs else ''))
 for label,top,bottom in [('single_vs_control','B7H3_CAR_T','T_cell'),('dual_vs_control','BC_CAR_T','T_cell'),('dual_vs_single','BC_CAR_T','B7H3_CAR_T')]:
  for d in ['B','E','M']:
   a=cpm.get(d+'_'+top);b=cpm.get(d+'_'+bottom);both=a is not None and b is not None and a>0 and b>0
   pairs.append(dict(gene=gene,comparison=label,reported_donor=d,numerator_CPM=a if a is not None else '',denominator_CPM=b if b is not None else '',log2_ratio=math.log2(a/b) if both else '',status='descriptive_both_positive' if both else 'zero_no_ratio' if counts else 'annotation_unresolved'))
 overview.append(dict(gene=gene,mapping_status=status,absent_published_CPM=gene in previous['genes_absent'],nonzero_samples=sum(v>0 for v in counts.values()) if counts else '',min_raw=min(counts.values()) if counts else '',max_raw=max(counts.values()) if counts else '',min_CPM=min(cpm.values()) if cpm else '',max_CPM=max(cpm.values()) if cpm else ''))
def write(n,rows):
 with (O/n).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for n,rows in [('gene_mapping.tsv',mapping),('all40_gene_sample_context.tsv',long),('condition_descriptive_summary.tsv',conditions),('within_donor_descriptive_ratios.tsv',pairs),('all40_overview.tsv',overview)]:write(n,rows)
assert len(long)==360 and len(conditions)==120 and len(pairs)==360
summary=dict(selected_genes=40,mapped=sum(x['mapping_status']=='exact_unique' for x in mapping),previously_missing_recovered=[x for x in overview if x['absent_published_CPM'] and x['mapping_status']=='exact_unique'],unresolved=[x for x in mapping if x['mapping_status']!='exact_unique'],all_zero=[x['gene'] for x in overview if x['mapping_status']=='exact_unique' and x['nonzero_samples']==0],input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'06_v1.22.7/gene_counts_unstranded.tsv',P/'06_v1.22.4/CAR_T_gene_expression_context.tsv']},limitations=['Reported donor labels, not genotype-verified identity','CPM denominator: all gene-assigned unstranded pairs','No p-values, no equivalence/stability inference, no pseudocount ratios','No change to candidate filters, rankings or genomic anchors','No expression equals accessibility or integration safety inference'])
(O/'context_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
