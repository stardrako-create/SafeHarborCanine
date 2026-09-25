"""Read-only validation and descriptive review after all nine corrected runs finish."""
from pathlib import Path
import csv,json,hashlib,statistics
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.7'
if not (O/'execution_complete.json').exists():
 print('Pending: corrected execution is not complete; no scientific results produced.');raise SystemExit(0)
done=json.loads((O/'execution_complete.json').read_text());assert done['samples']==9 and done['genes']==42309
manifest=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'));names=[r['column'] for r in manifest]
for name,h in done['matrix_sha256'].items():assert hashlib.sha256((O/name).read_bytes()).hexdigest()==h
individual={};qc=[]
for r in manifest:
 d=O/'samples'/r['SRR'];receipt=json.loads((d/'validated.json').read_text());assert receipt['signature']['sample']==r['column']
 for name,h in receipt['output_sha256'].items():assert hashlib.sha256((d/name).read_bytes()).hexdigest()==h
 lines=[l.split('\t') for l in (d/'STAR_ReadsPerGene.out.tab').read_text().splitlines()];individual[r['column']]={v[0]:list(map(int,v[1:])) for v in lines[4:]}
 totals=[sum(v[i] for v in individual[r['column']].values()) for i in range(3)];n=int(receipt['QC']['Number of input reads']);assert all(sum(int(v[i]) for v in lines)==n for i in [1,2,3])
 qc.append(dict(sample=r['column'],input_pairs=n,unique_pct=receipt['QC']['Uniquely mapped reads %'],assigned_unstranded=totals[0],assigned_fraction=totals[0]/n,forward_fraction_among_oriented=totals[1]/sum(totals[1:])))
for i,label in enumerate(['unstranded','forward','reverse']):
 with (O/f'gene_counts_{label}.tsv').open() as f:
  reader=csv.reader(f,delimiter='\t');assert next(reader)[1:]==names;seen=set()
  for v in reader:
   assert v[0] not in seen;seen.add(v[0]);assert list(map(int,v[1:]))==[individual[n][v[0]][i] for n in names]
  assert len(seen)==42309
genes=['ANO2','NTF3','NPNT','TBCK','TSEN15','C7H1orf21','CD247','CD8B','ARPC5','TET2']
context=[]
for g in genes:
 for row in qc:
  n=row['sample'];count=individual[n].get('gene-'+g,[None])[0]
  context.append(dict(gene=g,sample=n,unstranded_raw_count=count,CPM_gene_assigned_unstranded=None if count is None else count*1e6/row['assigned_unstranded']))
for filename,rows in [('review_sample_QC.tsv',qc),('review_focus_gene_context.tsv',context)]:
 with (O/filename).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
result=dict(status='independent_output_consistency_passed_scientific_limitations_remain',samples=9,genes=42309,comparison=json.loads((O/'comparison_with_v1226.json').read_text()),limitations=['CPM denominator is gene-assigned unstranded pairs, not necessarily the published CPM denominator','No differential expression or sample identity resolution claimed','Forward/reverse balance is suggestive, not complete library protocol validation','GEO declares SMART-Seq v4 3 prime DE; processing history needs reconciliation','RNA expression is not interval accessibility or integration safety'])
(O/'independent_review.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
