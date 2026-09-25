from pathlib import Path
import csv,gzip,json,math,statistics,hashlib,collections
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.4';old=P/'06_v1.22.3';samples={};current=None
with gzip.open(O/'GSE247355_family.soft.gz','rt') as f:
 for line in f:
  line=line.rstrip('\n')
  if line.startswith('^SAMPLE = '):current=line.split(' = ',1)[1];samples[current]=collections.defaultdict(list)
  elif line.startswith('^'):current=None
  elif current and line.startswith('!Sample_') and ' = ' in line:
   k,v=line.split(' = ',1);samples[current][k.removeprefix('!Sample_')].append(v)
matrixfile=old/'GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz'
with gzip.open(matrixfile,'rt') as f:reader=csv.DictReader(f,delimiter='\t');columns=reader.fieldnames;matrix=list(reader)
meta=[]
for acc,s in samples.items():
 ch=dict(x.split(': ',1) for x in s['characteristics_ch1'] if ': ' in x);title=s['title'][0];column=title.removesuffix('_1');assert column in columns[1:]
 r=dict(accession=acc,column=column,title=title,biological_replicate=ch.get('biological replicate'),treatment=ch.get('treatment'),cell_type=ch.get('cell type'),breed=ch.get('breed'),source_name=';'.join(s['source_name_ch1']),tissue=ch.get('tissue'),relations=';'.join(s['relation']),metadata_conflict='source_or_tissue_inconsistent_with_title_and_characteristics')
 expected_treatment='Untransduced' if column.endswith('_T_cell') else 'B7H3_CAR' if '_B7H3_' in column else 'BC_CAR'
 assert r['treatment']==expected_treatment
 specific=[v for v in [r['source_name'],r['tissue']] if v and v!='T cell']
 r['metadata_conflict']='specific_source_or_tissue_conflicts' if any(v!=column for v in specific) else 'specific_fields_agree' if specific else 'source_generic_tissue_missing'
 assert r['biological_replicate']==column.split('_')[0];meta.append(r)
assert len(meta)==9 and len({r['column'] for r in meta})==9
def write(name,rows):
 with (O/name).open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write('sample_metadata_audit.tsv',meta)
(O/'GEO_sample_metadata.json').write_text(json.dumps(samples,indent=2))
print(json.dumps(meta,indent=2))
lookup={r['Gene']:r for r in matrix};assert len(lookup)==14385
# Select genes from the existing annotated neighborhood, independently of their expression.
regions={'ANO2/NTF3':('NC_051831.1',39728324,39729324),'LOC119876429/LOC119872513':('NC_051811.1',17255706,17256706),'NPNT/TBCK':('NC_051836.1',27041161,27042161)}
genes=[]
for l in (P/'05_SHIP/canine_all_genes_stranded.bed').read_text().splitlines():
 v=l.split('\t')
 if len(v)>=6:genes.append(v)
context=[];selected=set(['ANO2','NTF3','LOC119876429','LOC119872513','NPNT','TBCK','TSEN15','C7H1orf21','ARPC5','TET2','CD3D','CD3E','CD3G','CD247','LCK','CD4','CD8A','CD8B'])
for name,(chrom,s,e) in regions.items():
 local=sorted([g for g in genes if g[0]==chrom],key=lambda g:max(int(g[1])-e,s-int(g[2]),0));seen=set()
 for g in local:
  if g[3] in seen:continue
  seen.add(g[3]);gap=max(int(g[1])-e,s-int(g[2]),0);context.append(dict(region=name,anchor_chrom=chrom,anchor_start=s,anchor_end=e,gene=g[3],gene_start=g[1],gene_end=g[2],strand=g[5],body_gap_bp=gap));selected.add(g[3])
  if len(seen)==10:break
write('local_gene_context.tsv',context)
summary=[];pairs=[];conditions={'T_control':'T_cell','B7H3_CAR':'B7H3_CAR_T','B7H3_CXCR2_CAR':'BC_CAR_T'}
for gene in sorted(selected):
 row=lookup.get(gene)
 for condition,suffix in conditions.items():
  values=[float(row[f'{d}_{suffix}']) for d in ['B','E','M']] if row else []
  summary.append(dict(gene=gene,condition=condition,in_CPM_table=row is not None,n_values=len(values),min_CPM=min(values) if values else '',median_CPM=statistics.median(values) if values else '',max_CPM=max(values) if values else '',nonzero_values=sum(x>0 for x in values) if values else '',B_CPM=values[0] if values else '',E_CPM=values[1] if values else '',M_CPM=values[2] if values else ''))
 for name,top,bottom in [('single_vs_control','B7H3_CAR_T','T_cell'),('dual_vs_control','BC_CAR_T','T_cell'),('dual_vs_single','BC_CAR_T','B7H3_CAR_T')]:
  for d in ['B','E','M']:
   a=float(row[f'{d}_{top}']) if row else None;b=float(row[f'{d}_{bottom}']) if row else None
   pairs.append(dict(gene=gene,comparison=name,biological_replicate=d,numerator_CPM=a,denominator_CPM=b,log2_ratio=math.log2(a/b) if a is not None and a>0 and b>0 else '',ratio_status='both_positive_descriptive_only' if a is not None and a>0 and b>0 else 'zero_value_no_pseudocount' if row else 'gene_absent_from_processed_table'))
write('CAR_T_gene_expression_context.tsv',summary);write('CAR_T_paired_descriptive_ratios.tsv',pairs)
validation=dict(status='passed',matrix_rows=len(matrix),matrix_columns=columns[1:],sample_metadata_matched=9,biological_replicate_labels=sorted({r['biological_replicate'] for r in meta}),selected_genes=len(selected),genes_present=sum(g in lookup for g in selected),genes_absent=[g for g in sorted(selected) if g not in lookup],gene_selection='Ten nearest distinct annotated gene symbols per declared region anchor plus predeclared candidate/marker genes; not expression-selected',methods='Median/range CPM across three reported replicate labels; within-label log2 ratios only when both values positive; no pseudocount or inferential testing',limitations=['Source/tissue fields in GEO may contradict specific sample title/treatment; titles/characteristics used and conflict preserved','No independent individual donor identity verification beyond biological replicate labels','Processed CPM is not raw count input for DESeq2','No local ATAC or post-integration safety inference','Transcript annotation and gene symbols differ between CanFam3.1 expression processing and ROS genomic context; missing symbols not zero'],sources={'GEO':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247355','paper':'https://pubmed.ncbi.nlm.nih.gov/38554158/','metadata_sha256':hashlib.sha256((O/'GSE247355_family.soft.gz').read_bytes()).hexdigest(),'CPM_sha256':hashlib.sha256(matrixfile.read_bytes()).hexdigest()})
(O/'RNA_context_validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation,indent=2));print(json.dumps([r for r in summary if r['gene'] in ['NPNT','TBCK','TSEN15','C7H1orf21','ARPC5','TET2']],indent=2))
