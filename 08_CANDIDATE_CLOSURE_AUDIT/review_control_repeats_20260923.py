from pathlib import Path
import os,json,csv,subprocess,datetime,hashlib,traceback
import pysam
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=C/'repeatmasker_context_20260923';O.mkdir(exist_ok=False)
def save(name,obj):(O/name).write_text(json.dumps(obj,indent=2))
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
status={'state':'running','started_utc':now()};save('status.json',status)
try:
 R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'
 rows=list(csv.DictReader((C/'transcript_and_boundary_review.tsv').open(),delimiter='\t'));contexts={}
 with pysam.FastaFile(str(R)) as ref,(O/'controls_context.fa').open('w') as f:
  for r in rows:
   s=int(r['start0']);e=int(r['end0']);lo=max(0,s-10000);hi=min(ref.get_reference_length(r['chrom']),e+10000)
   seq=ref.fetch(r['chrom'],lo,hi).upper();contexts[r['id']]={'chrom':r['chrom'],'context_start0':lo,'context_end0':hi,'window_start0':s,'window_end0':e,'sha256':hashlib.sha256(seq.encode()).hexdigest()};f.write('>'+r['id']+'\n'+seq+'\n')
 save('contexts.json',contexts)
 env=os.environ.copy();env['PATH']='/home/stardrako/miniforge3/envs/cactus/bin:'+env['PATH'];env['FAMDB_DATA_DIR']='/home/stardrako/miniforge3/envs/cactus/share/RepeatMasker/Libraries/famdb'
 cmd=['/home/stardrako/miniforge3/envs/cactus/bin/RepeatMasker','-species','dog','-pa','2','-dir',str(O),str(O/'controls_context.fa')];save('command.json',cmd)
 with (O/'stdout.log').open('w') as out,(O/'stderr.log').open('w') as err:p=subprocess.run(cmd,env=env,stdout=out,stderr=err)
 assert p.returncode==0,'RepeatMasker exit '+str(p.returncode)
 output=O/'controls_context.fa.out';assert output.exists()
 hits={k:[] for k in contexts}
 for line in output.read_text().splitlines():
  v=line.split()
  if len(v)<14 or not v[0].isdigit():continue
  key=v[4];c=contexts[key];start=c['context_start0']+int(v[5])-1;end=c['context_start0']+int(v[6]);lo=max(start,c['window_start0']);hi=min(end,c['window_end0'])
  if hi>lo:hits[key].append({'start0':lo,'end0':hi,'repeat':v[9],'class':v[10]})
 result=[]
 for key,c in contexts.items():
  mask=set()
  for h in hits[key]:mask.update(range(h['start0'],h['end0']))
  assert len(mask)<=1000
  result.append({'id':key,'masked_bp_union':len(mask),'fraction':len(mask)/1000,'hits':hits[key]})
 save('review.json',{'utc':now(),'results':result,'method':'RepeatMasker dog, default sensitivity; 10kb flanks then clip and union annotations inside original 1kb windows','output_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'limitations':['Repeat annotations depend on library and sensitivity; zero detected does not prove absence','No new candidate thresholds introduced']})
 status['state']='completed'
except Exception as e:
 status['state']='failed';status['error']=str(e);(O/'error.txt').write_text(traceback.format_exc())
finally:status['updated_utc']=now();save('status.json',status);print(json.dumps(status))
