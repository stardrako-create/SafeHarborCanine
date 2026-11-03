from pathlib import Path
import sys,subprocess,os,json,datetime,hashlib,shutil,traceback,collections
import pysam
A=Path(__file__).resolve().parent
O=A/'SV_PBGV000010/alignment_chunked/whole_genome_competition_20260926'
B='/home/stardrako/miniforge3/envs/atac/bin/bwa'
R=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa')
old=O.parent/'junction_competition_20260922_v2'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
if '--launch' in sys.argv:
 O.mkdir(exist_ok=False)
 with (O/'worker.log').open('w') as f:
  p=subprocess.Popen([sys.executable,str(Path(__file__).resolve())],stdout=f,stderr=f,start_new_session=True)
 print(json.dumps({'pid':p.pid,'output':str(O)}));sys.exit()
s={'state':'running','pid':os.getpid(),'started_utc':now()}
def save(stage):
 s.update(stage=stage,updated_utc=now());(O/'status.json').write_text(json.dumps(s,indent=2))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def run(cmd,out):
 with (O/out).open('w') as f,(O/(out+'.stderr.log')).open('w') as e:
  p=subprocess.Popen(cmd,stdout=f,stderr=e);s['child_pid']=p.pid;save(s['stage'])
  if p.wait()!=0:raise RuntimeError('Command failed: '+out)
 s['child_pid']=None
try:
 save('prepare_full_reference')
 refs=pysam.FastaFile(str(R));ref_lengths=dict(zip(refs.references,refs.lengths));refs.close()
 junctions=[(x.name,x.sequence) for x in pysam.FastxFile(str(old/'targets.fa')) if x.name.startswith('junction_')]
 assert len(junctions)==9 and not set(x[0] for x in junctions)&set(ref_lengths)
 combined=O/'genome_plus_junctions.fa'
 with combined.open('wb') as f:
  with R.open('rb') as src:shutil.copyfileobj(src,f,8*1024*1024)
  for n,seq in junctions:f.write(('\n>'+n+'\n'+seq+'\n').encode())
 pysam.faidx(str(combined))
 with pysam.FastaFile(str(combined)) as f:
  assert all(f.get_reference_length(n)==v for n,v in ref_lengths.items())
 reads={r.name:r.sequence for r in pysam.FastxFile(str(old/'reads.fa'))}
 assert len(reads)==30
 groups=sorted({n.rsplit('_mate',1)[0] for n in reads});assert len(groups)==15
 shutil.copyfile(old/'reads.fa',O/'reads.fa')
 for mate in (1,2):
  with (O/f'mate{mate}.fa').open('w') as f:
   for n in groups:f.write('>'+n+'\n'+reads[n+f'_mate{mate}']+'\n')
 (O/'inputs.json').write_text(json.dumps({'reference':str(R),'reference_sha256':sha(R),'combined_sha256':sha(combined),'reference_contigs':len(ref_lengths),'reference_bp':sum(ref_lengths.values()),'junctions':dict(junctions),'reads_sha256':sha(O/'reads.fa'),'pairs':15},indent=2))
 save('index_full_reference_plus_junctions');run([B,'index',str(combined)],'index.stdout')
 save('single_read_competition');run([B,'mem','-a','-t','4',str(combined),str(O/'reads.fa')],'single.sam')
 save('paired_read_competition');run([B,'mem','-a','-t','4',str(combined),str(O/'mate1.fa'),str(O/'mate2.fa')],'paired.sam')
 save('review')
 results={}
 for mode in ['single','paired']:
  by=collections.defaultdict(list)
  with pysam.AlignmentFile(str(O/(mode+'.sam'))) as f:
   for r in f:
    if r.is_unmapped:continue
    j=r.reference_name.startswith('junction_');blocks=r.get_blocks()
    left=sum(max(0,min(e,1000)-b) for b,e in blocks if b<1000)
    right=sum(max(0,e-max(b,1000)) for b,e in blocks if e>1000)
    name=r.query_name+(('_mate1' if r.is_read1 else '_mate2') if mode=='paired' else '')
    by[name].append(dict(target=r.reference_name,start0=r.reference_start,end0=r.reference_end,AS=r.get_tag('AS') if r.has_tag('AS') else None,mapq=r.mapping_quality,cigar=r.cigarstring,reverse=r.is_reverse,secondary=r.is_secondary,supplementary=r.is_supplementary,proper_pair=r.is_proper_pair,mate_reference=r.next_reference_name,template_length=r.template_length,junction=j,spanning=j and min(left,right)>=20))
  summaries=[]
  for n in reads:
   hits=by[n];ja=max((h['AS'] for h in hits if h['spanning'] and h['AS'] is not None),default=None);ga=max((h['AS'] for h in hits if not h['junction'] and h['AS'] is not None),default=None)
   summaries.append(dict(read=n,best_spanning_AS=ja,best_genomic_AS=ga,margin=None if ja is None or ga is None else ja-ga))
  (O/(mode+'_hits.json')).write_text(json.dumps(by,indent=2));results[mode]=summaries
 report={'utc':now(),'results':results,'limitations':['Selected 15 pairs, not unbiased independent SV evidence','BWA heuristic reported alignments, not exhaustive enumeration','Paired mode uses default BWA estimation on only 15 selected pairs; inspect stderr for insufficient estimation; no calibrated pair likelihood','Augmented-reference MAPQ is not original-reference MAPQ','AS differences do not establish genotype or resolve repetitive origin','Same source reads as prior panel; no increase in independent evidence']}
 (O/'review.json').write_text(json.dumps(report,indent=2));s.update(state='completed',review_sha256=sha(O/'review.json'));save('completed')
except Exception:
 s['state']='failed';(O/'error.txt').write_text(traceback.format_exc());save('failed');raise
