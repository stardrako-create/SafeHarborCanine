from pathlib import Path
import json,numpy as np,pyBigWig,subprocess,shutil,time,hashlib
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.8';D=O/'ATAC'
manifest=json.loads((D/'build_manifest.json').read_text());h=pyBigWig.open(str(D/'mother_track_accessibility_level.bw'));chroms=h.chroms();h.close();sums=dict(mean=0,var=0,evidence=0);bins=0;vals=[]
for c,size in chroms.items():
 with np.load(D/'cache'/(c+'.npz')) as d:
  m=d['mean'];v=d['variability'];n=d['evidence'];assert len(m)==(size+24)//25 and len(v)==len(n)==len(m)
  assert np.isfinite(m).all() and (m>=0).all() and (n<=76).all();assert np.array_equal(np.isfinite(v),n>=2);assert (v[n>=2]>=0).all();assert (m[n==0]==0).all();assert (m[n>0]>0).all()
  lengths=np.minimum(25,size-np.arange(len(m))*25);sums['mean']+=int(lengths[n>0].sum());sums['var']+=int(lengths[n>=2].sum());sums['evidence']+=size;bins+=len(n);vals.append(m[n>0])
median=float(np.median(np.concatenate(vals)));assert median==manifest['gain_background'];del vals
for k,f in [('mean','mother_track_accessibility_level.bw'),('var','variability.bw'),('evidence','evidence_dogs.bw')]:
 with pyBigWig.open(str(D/f)) as h:assert h.header()['nBasesCovered']==sums[k],(k,h.header(),sums[k])
(O/'rebuild_global_validation.json').write_text(json.dumps(dict(status='passed',sequences=len(chroms),bins_checked=bins,covered_bp=sums,reproduced_gain_background=median,scope='All cached bins and final BigWig header coverage; independent base-level numerical checks recorded separately.'),indent=2));print('Global cache validation passed',bins,flush=True)
code=O/'code';code.mkdir(exist_ok=True)
for n in ['score_ship_candidates_v2.py','bw_utils.py','build_mother_track_v2.py']:shutil.copy2(P/'scripts'/n,code/n)
cmd=json.loads((P/'06_v1.21.5/scoring_command.json').read_text());cmd[1]=str(code/'score_ship_candidates_v2.py')
for flag,path in [('--atac-mean-bw',D/'mother_track_accessibility_level.bw'),('--atac-variability-bw',D/'variability.bw'),('--out-scored',O/'candidates_scored_v1218.tsv'),('--out-passing-bed',O/'passing_legacy_checks_v1218.bed')]:cmd[cmd.index(flag)+1]=str(path)
(O/'scoring_command.json').write_text(json.dumps(cmd,indent=2));t=time.time()
with (O/'scoring.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
(O/'scoring_execution.json').write_text(json.dumps(dict(returncode=r.returncode,seconds=time.time()-t),indent=2));print('Scoring exit',r.returncode,flush=True);assert r.returncode==0
