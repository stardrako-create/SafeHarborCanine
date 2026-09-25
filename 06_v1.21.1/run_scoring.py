from pathlib import Path
import csv,json,subprocess,sys,hashlib,collections,shutil,datetime
import pysam,pyBigWig,numpy as np
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'06_v1.21.1'
def read(f):
 with f.open(encoding='utf-8') as x:return list(csv.DictReader(x,delimiter='\t'))
def write(f,rows):
 with f.open('w',newline='',encoding='utf-8') as x:
  w=csv.DictWriter(x,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def key(r):return r['chrom'],int(r['start']),int(r['end'])
manifest=json.loads((P/'06_v1.21.0/manifest.json').read_text())
rows=read(P/'05_SHIP/ship_raw_candidates.tsv')
fa_path=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'
fa=pysam.FastaFile(str(fa_path))
repeat=[]
for r in rows:
 c,s,e=key(r);seq=fa.fetch(c,s,e)
 masked=sum(ch in 'acgt' for ch in seq)
 repeat.append({'chrom':c,'start':s,'end':e,'length':e-s,'masked_bp':masked,'pct_repeat':round(100*masked/(e-s),2)})
fa.close()
prior={key(r):r for r in read(P/'05_SHIP/repeat_content_v5candidates.tsv')}
comp=[{'chrom':r['chrom'],'start':r['start'],'end':r['end'],'old':float(prior[key(r)]['pct_repeat']),'softmask':r['pct_repeat']} for r in repeat if key(r) in prior]
write(O/'repeat_softmask_validation.tsv',comp)
maxdiff=max(abs(r['old']-r['softmask']) for r in comp)
print('Repeat softmask maximum difference from existing annotations:',maxdiff,flush=True)
write(O/'repeat_softmask_all_candidates.tsv',repeat)
# Only promote soft-mask counts when they exactly reproduce the established metric.
if maxdiff<=0.01:
 manifest['inputs']['repeat-content-tsv']['path']=str(O/'repeat_softmask_all_candidates.tsv')
else:
 print('Softmask differs: keeping prior repeat annotations, new measurements diagnostic only',flush=True)
cmd=[sys.executable,str(O/'code/score_ship_candidates_v2.py')]
for arg,info in manifest['inputs'].items():cmd+=['--'+arg,info['path']]
cmd+=['--min-atac-accessibility-percentile','0.55','--out-scored',str(O/'candidates_scored_v1211.tsv'),'--out-passing-bed',str(O/'candidates_passing_recorded_checks.bed')]
with (O/'scoring.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
scored=read(O/'candidates_scored_v1211.tsv');old={key(r):r for r in read(P/'06_v1.20.0/candidates_scored_v120d.tsv')}
assert len(scored)==461
eligible=sorted([r for r in scored if r['hard_veto']=='False'],key=lambda r:float(r['final_score']) if r['final_score'] else -1,reverse=True)
write(O/'provisional_candidates.tsv',eligible)
comparison=[]
for r in eligible:
 prev=old[key(r)]
 comparison.append({'chrom':r['chrom'],'start':r['start'],'end':r['end'],'genes':r['left_gene']+'/'+r['right_gene'],'previous_score':prev['final_score'],'corrected_score':r['final_score'],'previous_rrbs':prev['rrbs_mean'],'corrected_rrbs':r['rrbs_mean'],'rrbs_observed_fraction':r['rrbs_observed_fraction'],'repeat_pct':r['pct_repeat'],'status':r['evaluation_status'],'missing':r['missing_evidence']})
write(O/'comparison.tsv',comparison)
assert all(r['evaluation_status']!='passes_recorded_checks' for r in scored if r['missing_evidence'])
for r in eligible:
 assert 0<float(r['rrbs_observed_fraction'])<=1
summary={'counts':dict(collections.Counter(r['evaluation_status'] for r in scored)),'repeat_max_difference':maxdiff,'comparison':comparison}
(O/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
run={'version':'1.21.1','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd,'scope':'corrected scoring over existing ATAC/RRBS tracks; RRBS mean masked by measured coverage for candidates and backgrounds; no ATAC rebuild','python':sys.version,'pyBigWig':pyBigWig.__version__,'inputs':{}}
for arg,info in manifest['inputs'].items():
 f=Path(info['path']);run['inputs'][arg]={'path':str(f),'size':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns}
 h=hashlib.sha256()
 with f.open('rb') as fin:
  for block in iter(lambda:fin.read(8*1024*1024),b''):h.update(block)
 run['inputs'][arg]['sha256']=h.hexdigest()
run['code_sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in (O/'code').glob('*.py')}
(O/'manifest.json').write_text(json.dumps(run,indent=2),encoding='utf-8')
shutil.copy2(Path(__file__),O/'run_scoring.py')
print(json.dumps(summary,indent=2),flush=True)
