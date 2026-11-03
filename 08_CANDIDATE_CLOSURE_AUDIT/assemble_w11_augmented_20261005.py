from pathlib import Path
import pysam,shutil,subprocess,os,json,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_AUGMENTED_ASSEMBLY_20261005';O.mkdir(exist_ok=False);W=Path('/mnt/d/SafeHarborWork/w11_augmented_20261005');W.mkdir(exist_ok=False);base=A/'W11_SPADES_20261004';raw=A/'W11_RAW_RECRUIT_20261005';names={'SRR14750332.53211536','SRR14750332.171739498'}
def state(s,**kw):(O/'status.json').write_text(json.dumps({'stage':s,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw},indent=2))
try:
 state('preparing')
 for m in [1,2]:
  p=W/f'R{m}.fastq';shutil.copy2(base/f'R{m}.fastq',p);n=0
  with p.open('a') as out,pysam.FastxFile(str(raw/f'selected_R{m}.fastq')) as f:
   for r in f:
    if r.name in names:out.write('@'+r.name+'/'+str(m)+'\n'+r.sequence+'\n+\n'+r.quality+'\n');n+=1
  assert n==2
 shutil.copy2(base/'single.fastq',W/'single.fastq')
 (O/'inputs.json').write_text(json.dumps({'pairs':7358,'singles':4,'added_names':sorted(names),'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in W.glob('*.fastq')},'limits':['Exploratory inclusion; new pairs not added to validated support counts']},indent=2))
 cmd=['python3','/mnt/d/SafeHarborWork/w11_spades_20261004_ascii/SPAdes-4.3.0-Linux/bin/spades.py','--only-assembler','-1',str(W/'R1.fastq'),'-2',str(W/'R2.fastq'),'-s',str(W/'single.fastq'),'-k','21,33,55,77','-t','6','-m','12','-o',str(W/'assembly')];env=os.environ.copy();env['LD_LIBRARY_PATH']='/home/stardrako/miniforge3/envs/atac/lib';(O/'command.json').write_text(json.dumps(cmd));state('assembling')
 with (W/'run.log').open('w') as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
 with pysam.FastaFile(str(A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa')) as f:j=f.fetch('junction_L4707939_R33613308').upper()
 with pysam.FastaFile(str(A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa')) as f:tail='A'*20+rc(f.fetch('NC_051836.1',35342033,35342133).upper())[:20]
 rows=[]
 with pysam.FastxFile(str(W/'assembly/contigs.fasta')) as f:
  for r in f:
   seq=r.sequence.upper();rows.append({'name':r.name,'length':len(seq),'junction':j[980:1020] in seq or j[980:1020] in rc(seq),'tail':tail in seq or tail in rc(seq)})
 (O/'motif_review.json').write_text(json.dumps(rows,indent=2));state('completed_needs_review',contigs=len(rows),motif_hits=[r for r in rows if r['junction'] or r['tail']])
except Exception as e:state('failed',error=str(e));raise