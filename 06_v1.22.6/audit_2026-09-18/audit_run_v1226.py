from pathlib import Path
import json,csv,gzip,collections,hashlib,shutil,re
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino'); O=P/'06_v1.22.6'; A=O/'audit_2026-09-18'; A.mkdir(exist_ok=True)
manifest=json.loads((O/'ENA_verified_download_manifest.json').read_text()); runs=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'))
checks=[]; quality=[]
for f in manifest['files']:
 p=O/'raw'/Path(f['path']).name; cert=json.loads(p.with_suffix(p.suffix+'.verified.json').read_text()); assert p.stat().st_size==f['bytes']==cert['bytes']; assert cert['observed_md5']==cert['expected_md5']==f['md5']
 checks.append(dict(file=p.name,bytes=p.stat().st_size,certificate_verified=True,md5_rehashed_this_audit=False))
 q=collections.Counter()
 with gzip.open(p,'rt') as g:
  for i in range(10000):
   head,seq,plus,qual=[g.readline().rstrip() for _ in range(4)]; assert head.startswith('@') and plus.startswith('+') and len(seq)==len(qual); q.update(ord(c)-33 for c in qual)
 quality.append(dict(file=p.name,first_reads=10000,phred_histogram=dict(q)))
indexids=[l.split('\t')[0] for l in (O/'reference/STAR_index/geneInfo.tab').read_text().splitlines()[1:]]
annot=json.loads((O/'reference/annotation_validation.json').read_text()); missing=sorted(set(annot['gene_ids'])-set(indexids)); (A/'annotation_genes_missing_from_STAR.json').write_text(json.dumps(missing,indent=2))
rows=[]; countsets={}; hashes={}
for r in runs:
 d=O/'samples'/r['SRR']; assert 'ALL DONE!' in (d/'STAR_Log.out').read_text(); assert 'finished successfully' in (d/'STAR_console.log').read_text()
 fp=json.loads((d/'fastp.json').read_text()); qc={k.strip():v.strip() for l in (d/'STAR_Log.final.out').read_text().splitlines() if '|' in l for k,v in [l.split('|',1)]}
 n=int(qc['Number of input reads']); assert fp['summary']['after_filtering']['total_reads']==2*n
 lines=[l.split('\t') for l in (d/'STAR_ReadsPerGene.out.tab').read_text().splitlines()]; assert [v[0] for v in lines[:4]]==['N_unmapped','N_multimapping','N_noFeature','N_ambiguous']; assert [v[0] for v in lines[4:]]==indexids
 assert all(len(v)==4 and all(x.isdigit() for x in v[1:]) for v in lines)
 assert all(sum(int(v[i]) for v in lines)==n for i in [1,2,3])
 totals=[sum(int(v[i]) for v in lines[4:]) for i in [1,2,3]]; countsets[r['column']]={v[0]:list(map(int,v[1:])) for v in lines[4:]}
 rows.append(dict(sample=r['column'],run=r['SRR'],input_pairs=n,unique_pct=float(qc['Uniquely mapped reads %'].strip('%')),assigned_unstranded=totals[0],assigned_pct=100*totals[0]/n,forward_share=totals[1]/sum(totals[1:]),removed_reads=fp['summary']['before_filtering']['total_reads']-2*n,q30=fp['summary']['before_filtering']['q30_rate']))
for i,label in enumerate(['unstranded','forward','reverse']):
 p=O/f'gene_counts_{label}.tsv'; reader=csv.reader(p.open(),delimiter='\t'); names=next(reader)[1:]; assert names==[r['column'] for r in runs]
 ng=0
 for v in reader:
  assert list(map(int,v[1:]))==[countsets[n][v[0]][i] for n in names]; ng+=1
 assert ng==len(indexids); hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
with (A/'sample_QC.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader();w.writerows(rows)
focus=['ANO2','NTF3','CD247','CD8B','NPNT','TBCK','TSEN15','C7H1orf21','ARPC5','TET2']
context={g:{n:countsets[n].get('gene-'+g) for n in countsets} for g in focus}
result=dict(status='execution_and_matrix_consistency_verified_scientific_review_open',samples=rows,missing_annotation_genes=missing,quality_sample=quality,input_certificate_checks=checks,matrix_sha256=hashes,focus_gene_raw_counts_all_orientations=context,disk_free_bytes=shutil.disk_usage(O).free)
(A/'audit.json').write_text(json.dumps(result,indent=2)); print(json.dumps(dict(samples=rows,missing_count=len(missing),missing_first=missing[:12],quality=[x['phred_histogram'] for x in quality],focus=context),indent=2))
