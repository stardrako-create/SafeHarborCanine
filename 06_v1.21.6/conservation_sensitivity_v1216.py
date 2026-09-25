from pathlib import Path
import subprocess,json,time
O=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/06_v1.21.6');Z=Path('/mnt/d/Jin2024_work/zoonomia_hal');B='/home/stardrako/miniforge3/envs/cactus/bin/'
cmd=[B+'phyloFit','--init-model',str(Z/'neutral_model_1region.mod'),'--scale-only','--no-freqs','--no-rates','--precision','MED','--msa-format','MAF','--out-root',str(O/'independent_region_scale'),str(Z/'maf_dog/neutral_2.chr1.97942942.maf')]
t=time.time()
with (O/'conservation_scale.log').open('w') as out:
 try:r=subprocess.run(cmd,stdout=out,stderr=subprocess.STDOUT,timeout=600);status=r.returncode
 except subprocess.TimeoutExpired:status='timeout600s'
(O/'conservation_scale_execution.json').write_text(json.dumps(dict(command=cmd,status=status,seconds=time.time()-t,limitation='Sensitivity check on a second random unfiltered region; branch proportions, frequencies and rates inherited from pilot. Not an independently validated neutral model.'),indent=2));print(status,flush=True)
