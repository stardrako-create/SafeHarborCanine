import pathlib,json,subprocess,hashlib,shutil,fcntl,datetime
A=pathlib.Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/08_CANDIDATE_CLOSURE_AUDIT')
P=A/'SV_PBGV000012';O=P/'ROS_integration_20261001';O.mkdir(exist_ok=True)
l=open(O/'worker.lock','w');fcntl.flock(l,fcntl.LOCK_EX|fcntl.LOCK_NB)
S='/home/stardrako/miniforge3/envs/atac/bin/samtools'
def state(x,**kw):(O/'status.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':x,**kw},indent=2))
def run(args,name):
 with open(O/(name+'.out'),'w') as out,open(O/(name+'.err'),'w') as err:subprocess.run([S,*map(str,args)],stdout=out,stderr=err,check=True)
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as inp:
  for b in iter(lambda:inp.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
try:
 state('auditing_chunks')
 files=[P/'ROS_chunked_v2/chunk_000001_1000000/receipt.json']+list((P/'ROS_production_20260930').glob('chunk_*/receipt.json'))
 rows=sorted([(json.loads(f.read_text()),f) for f in files],key=lambda x:x[0]['start']);end=0;selected=0;manifest=[]
 assert len(rows)==215
 for i,(r,f) in enumerate(rows):
  assert r['start']==end+1 and r['pairs']==r['end']-r['start']+1 and r['primary_records']==2*r['pairs'];end=r['end']
  b=f.parent/'selected.bam';idx=f.parent/'selected.bam.bai'
  h=sha(b);hi=sha(idx)
  if i>0:
   saved=json.loads((f.parent/'sha256.json').read_text());assert saved['selected.bam']==h and saved['selected.bam.bai']==hi
  run(['quickcheck','-v',b],f'quickcheck_{i}')
  count=int(subprocess.check_output([S,'view','-c',str(b)]));assert count==r['selected_records'];selected+=count
  manifest.append({'bam':str(b),'sha256':h,'index_sha256':hi,'count':count,'start':r['start'],'end':end})
 assert end==214126885
 (O/'input_audit.json').write_text(json.dumps({'pairs':end,'selected_records':selected,'chunks':manifest,'limitation':'selected reference-or-mate chr32 alignments only; not whole-genome BAM'},indent=2))
 needed=sum(pathlib.Path(x['bam']).stat().st_size for x in manifest)*8+5*1024**3
 assert shutil.disk_usage(O).free>needed,'Insufficient integration disk budget'
 assert not (O/'merged.bam').exists(),'Preserve existing integration outputs'
 (O/'inputs.list').write_text(''.join(x['bam']+'\n' for x in manifest))
 state('merging',pairs=end,selected_records=selected)
 run(['merge','-@','2','-c','-p','-b',O/'inputs.list',O/'merged.bam'],'merge')
 run(['quickcheck','-v',O/'merged.bam'],'merged_check')
 assert int(subprocess.check_output([S,'view','-c',str(O/'merged.bam')]))==selected
 state('name_sort_fixmate')
 run(['sort','-n','-@','2','-m','512M','-o',O/'namesort.bam',O/'merged.bam'],'namesort')
 run(['fixmate','-m',O/'namesort.bam',O/'fixmate.bam'],'fixmate')
 run(['sort','-@','2','-m','512M','-o',O/'coordinate.bam',O/'fixmate.bam'],'coordinate')
 state('marking_duplicates')
 run(['markdup','-s',O/'coordinate.bam',O/'selected_markdup.bam'],'markdup')
 run(['quickcheck','-v',O/'selected_markdup.bam'],'final_check')
 run(['index',O/'selected_markdup.bam'],'index')
 run(['flagstat',O/'selected_markdup.bam'],'flagstat')
 count=int(subprocess.check_output([S,'view','-c',str(O/'selected_markdup.bam')]))
 assert count==selected
 (O/'completion.json').write_text(json.dumps({'pairs':end,'chunks':len(rows),'alignment_records':count,'bam_sha256':sha(O/'selected_markdup.bam'),'duplicates':'marked globally, not removed','scope':'chr32 linked subset; structural review pending'},indent=2))
 state('completed_structural_review_pending')
except Exception as e:state('failed',error=str(e));raise