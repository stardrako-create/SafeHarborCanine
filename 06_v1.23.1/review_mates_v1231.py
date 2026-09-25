from pathlib import Path
import json,csv,hashlib,numpy as np
from scipy.stats import spearmanr
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.23.1'
done=json.loads((O/'execution_complete.json').read_text());assert done['alignments']==27 and done['samples']==9
qc=list(csv.DictReader((O/'mate_alignment_QC.tsv').open(),delimiter='\t'));corr=list(csv.DictReader((O/'mate_count_concordance.tsv').open(),delimiter='\t'));assert len(qc)==len(corr)==27
ids=[l.split('\t')[0] for l in (P/'06_v1.22.7/STAR_index/geneInfo.tab').read_text().splitlines()[1:]];vectors={};checked=[]
for row in qc:
 d=O/row['run']/row['mode'];r=json.loads((d/'validated.json').read_text());p=d/'STAR_ReadsPerGene.out.tab';assert hashlib.sha256(p.read_bytes()).hexdigest()==r['counts_sha256']
 assert 'ALL DONE!' in (d/'STAR_Log.out').read_text() and 'finished successfully' in (d/'console.log').read_text()
 for i,h in enumerate(r['signature']['input_sha256'],1):assert hashlib.sha256((O/row['run']/f'first100k_R{i}.fastq').read_bytes()).hexdigest()==h
 data=[l.split('\t') for l in p.read_text().splitlines()];assert [v[0] for v in data[4:]]==ids and len(ids)==42309
 assert all(sum(int(v[i]) for v in data)==100000 for i in [1,2,3]);v=np.array([int(x[1]) for x in data[4:]]);vectors[(row['sample'],row['mode'])]=v
 assert int(v.sum())==int(row['assigned_unstranded']) and int((v>0).sum())==int(row['detected_genes'])
 assert abs(float(row['assigned_pct'])-v.sum()/1000)<1e-10
 actual={k.strip():val.strip() for l in (d/'STAR_Log.final.out').read_text().splitlines() if '|' in l for k,val in [l.split('|',1)]};assert actual==r['QC'] and float(actual['Uniquely mapped reads %'].strip('%'))==float(row['unique_pct'])
 checked.append(row['run']+'/'+row['mode'])
for row in corr:
 a=vectors[(row['sample'],row['mode_a'])];b=vectors[(row['sample'],row['mode_b'])];m=(a+b)>0;assert int(m.sum())==int(row['genes_nonzero_either']);assert abs(spearmanr(a[m],b[m]).statistic-float(row['spearman_nonzero_either']))<1e-12
summary={'status':'independent_review_passed','validated_alignments':checked,'mode_ranges':{mode:{key:[min(float(r[key]) for r in qc if r['mode']==mode),max(float(r[key]) for r in qc if r['mode']==mode)] for key in ['unique_pct','assigned_pct']} for mode in ['paired','R1','R2']},'correlation_ranges':{a+'_'+b:[min(float(r['spearman_nonzero_either']) for r in corr if r['mode_a']==a and r['mode_b']==b),max(float(r['spearman_nonzero_either']) for r in corr if r['mode_a']==a and r['mode_b']==b)] for a,b in [('R1','paired'),('R2','paired'),('R1','R2')]}}
(O/'independent_review.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
