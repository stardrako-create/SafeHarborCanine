from pathlib import Path
import csv,json,hashlib,shutil,collections
import numpy as np,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.5'
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
old=read(P/'06_v1.21.3/candidates_scored_v1213.tsv');new=read(O/'candidates_scored_coherent76.tsv');assert len(old)==len(new)==461
key=lambda r:(r['chrom'],r['start'],r['end']);a={key(r):r for r in old};b={key(r):r for r in new};assert a.keys()==b.keys()
for k in a:
 for field in ['final_score','veto_atac_peak','evaluation_status']:assert a[k][field]==b[k][field],(k,field)
per=read(O/'per_dog_local_atac.tsv');assert len(per)==456
for gene in {r['genes'] for r in per}:
 rs=[r for r in per if r['genes']==gene];assert len({r['sample'] for r in rs})==76
 for r in rs:assert 0<=float(r['gate_fraction'])<=1 and float(r['cpm_mean'])>=0
tiles=read(O/'local_remap_tiles.tsv');assert len(tiles)==309;passed=sum(r['passes_remap_proxy']=='True' for r in tiles)
epic=read(O/'epic_local_states.tsv');coverage=collections.Counter()
for r in epic:coverage[r['genes'],r['tissue_code']]+=int(r['observed_bp']);assert 1<=int(r['state'])<=13
assert len(coverage)==33 and all(0<n<=1000 for n in coverage.values())
vcf=pysam.VariantFile('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz')
v=next(v for v in vcf.fetch('chr32',13062402,13063402) if v.id=='chr32:4707939:DG')
carriers=[s for s in v.samples.values() if any(x not in (0,None) for x in (s.get('GT') or []))]
sv={'id':v.id,'start':v.start,'end':v.stop,'span_bp':v.stop-v.start,'AC':v.info['AC'],'AN':v.info['AN'],'AF':v.info['AF'],'carrier_genotypes':len(carriers),'PASS_carriers':sum(s.get('FT')=='PASS' for s in carriers),'interpretation':'Rare catalogued large deletion overlaps NPNT local window; not verified in project animals, no automatic veto. SNPs/indels not assessed.'}
(O/'NPNT_SV_review.json').write_text(json.dumps(sv,indent=2))
qc=read(O/'qc_by_dog.tsv')
stats={'validation_date':'2026-09-14','coherent76_score_rows':461,'score_status_and_atac_veto_identical_to_v1213':True,'statuses':dict(collections.Counter(r['evaluation_status'] for r in new)),'per_dog_measurements':len(per),'samples':76,'local_read_errors':0,'remap_tiles':309,'remap_pass':passed,'epic_tissue_window_pairs':33,'epic_states_observed':sorted({int(r['state']) for r in epic}),'epic_coverage_by_window_tissue':{g+'|'+t:n for (g,t),n in coverage.items()},'qc_frip_at_least_0_2':sum(float(r['frip'])>=.2 for r in qc),'qc_note':'FRiP computed against per-sample called peaks; custom TSS metric averages broad windows, not equivalent to ENCODE TSS-peak enrichment. No absolute QC pass claimed.','NPNT_SV':sv}
(O/'validation.json').write_text(json.dumps(stats,indent=2))
for filename in ['annotations_v1215.py','validate_clean_v1215.py']:
 shutil.copy2(Path('/mnt/c/Users/Utilizador/OneDrive/Documents/New project')/filename,O/filename)
(O/'file_hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='file_hashes.json'},indent=2))
print(json.dumps(stats,indent=2))
