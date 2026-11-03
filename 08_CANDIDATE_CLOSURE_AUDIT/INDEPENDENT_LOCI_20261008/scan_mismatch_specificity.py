from pathlib import Path
import json,csv,pysam,subprocess,hashlib,os,datetime,gzip,collections,random,fcntl,traceback
N=Path(__file__).resolve().parent;P=N.parent.parent;O=N/'mismatch_specificity_20261009';O.mkdir(exist_ok=False)
lock=open('/tmp/safeharbor_mismatch_scan_20261009.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
B=Path('/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin');R=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';idx=Path('/mnt/d/Jin2024_work/bowtie1_ROS_20261009');idx.mkdir(exist_ok=False)
status=dict(state='preparing',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reference=str(R),index=str(idx/'ROS'),threads=8,alignments_checked=0)
def save():
 temp=O/'status.tmp'
 with temp.open('w') as f:json.dump(status,f,indent=2);f.flush();os.fsync(f.fileno())
 temp.replace(O/'status.json')
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for part in iter(lambda:f.read(8388608),b''):h.update(part)
 return h.hexdigest()
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def run(cmd,stem,seconds):
 (O/(stem+'.command.json')).write_text(json.dumps(cmd))
 with (O/(stem+'.out')).open('w') as out,(O/(stem+'.err')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=seconds)
try:
 save();rng=random.Random(79);toy=''.join(rng.choice('ACGT') for _ in range(160));toy+='A'*30+toy[10:40];toyref=O/'toy.fa';toyref.write_text('>toy\n'+toy+'\n');q={}
 for i,s in enumerate([toy[:20],toy[20:40],rc(toy[40:60]),'A'*20]):
  q[str(i)]=s
  for k in [1,2,3]:
   t=list(s)
   for z in range(k):t[z]='A' if t[z]!='A' else 'C'
   q[f'{i}_{k}']=''.join(t)
 toyqueries=O/'toy_queries.fa';toyqueries.write_text(''.join(f'>{k}\n{v}\n' for k,v in q.items()));run([str(B/'bowtie-build'),'--quiet',str(toyref),str(O/'toy')],'toy_build',120)
 run([str(B/'bowtie'),'-f','-v','3','-a','--quiet',str(O/'toy'),str(toyqueries)],'toy_scan',120)
 actual=set()
 for line in (O/'toy_scan.out').read_text().splitlines():
  x=line.split('\t');actual.add((x[0],x[1],int(x[3])))
 expected={(key,strand,i) for key,seq in q.items() for strand,query in [('+',seq),('-',rc(seq))] for i in range(len(toy)-19) if sum(a!=b for a,b in zip(query,toy[i:i+20]))<=3}
 assert actual==expected,('Bowtie -a -v3 differs from brute force',len(actual),len(expected));status['self_test']=dict(queries=len(q),alignments=len(actual),matches_brute_force=True)
 status['reference_sha256']=digest(R);status['tool_sha256']={f.name:digest(f) for f in [B/'bowtie-align-s',B/'bowtie-build-s']};status['state']='building_index';save()
 run([str(B/'bowtie-build'),'--threads','8',str(R),str(idx/'ROS')],'index_build',7200)
 indexfiles=sorted(idx.glob('ROS.*.ebwt'));assert len(indexfiles)==6;status['index_files']=[dict(file=str(f),size=f.stat().st_size,sha256=digest(f)) for f in indexfiles]
 masks={r['window']:set(r['mask_offsets']) for r in json.loads((N/'two_loci_REF_mask_specificity.json').read_text())['results']};A=N.parent/'ALL_CANDIDATES_20261005';ano='NC_051831.1_39708074_39709074';masks[ano]=set(next(r for r in json.loads((A/'masked_specificity_20261008.json').read_text())['results'] if r['id']==ano)['masked_offsets']);fa=pysam.FastaFile(str(R));queries={};sourcefiles=[N/'two_loci_REF_mask_specificity.json',A/'masked_specificity_20261008.json']
 for tsv,field in [(N/'exact_specificity_next9/all_short_segments_exact_uniqueness.tsv','window_id'),(A/'exact_specificity_20261006/all_short_segments_exact_uniqueness.tsv','window_id')]:
  sourcefiles.append(tsv)
  with tsv.open() as f:
   for r in csv.DictReader(f,delimiter='\t'):
    key=r[field]
    if key not in masks or int(r['length'])!=20 or r['unique_exact']!='True':continue
    chrom,start,end=key.rsplit('_',2);lo=int(r['start0']);off=lo-int(start)
    if any(i in masks[key] for i in range(off,off+20)):continue
    name=key+'__'+str(off);queries[name]=dict(window=key,chrom=chrom,start0=lo,end0=lo+20,offset0=off,sequence=fa.fetch(chrom,lo,lo+20).upper(),self_hits=0,off_target_counts={str(i):0 for i in range(4)},first10_offtargets=[])
 assert len(queries)==2127,len(queries)
 status['source_hashes']={str(f):digest(f) for f in sourcefiles};queryfile=O/'queries.fa';queryfile.write_text(''.join(f'>{k}\n{r["sequence"]}\n' for k,r in queries.items()));status.update(state='scanning',queries=len(queries));save()
 cmd=[str(B/'bowtie'),'-f','-v','3','-a','--threads','4',str(idx/'ROS'),str(queryfile)];(O/'scan.command.json').write_text(json.dumps(cmd));seen=set()
 with (O/'scan.err').open('w') as err,gzip.open(O/'alignments.tsv.gz','wt',compresslevel=1) as raw:
  proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True)
  for line in proc.stdout:
   raw.write(line);x=line.rstrip('\n').split('\t');assert len(x)>=7;xkey=x[0];r=queries[xkey];strand=x[1];chrom=x[2];lo=int(x[3]);matchkey=(xkey,strand,chrom,lo);assert matchkey not in seen;seen.add(matchkey)
   target=fa.fetch(chrom,lo,lo+20).upper();oriented=r['sequence'] if strand=='+' else rc(r['sequence']);nm=sum(a!=b for a,b in zip(oriented,target));assert len(target)==20 and nm<=3
   if chrom==r['chrom'] and lo==r['start0']:assert nm==0;r['self_hits']+=1
   else:
    r['off_target_counts'][str(nm)]+=1
    if len(r['first10_offtargets'])<10:r['first10_offtargets'].append(dict(chrom=chrom,start0=lo,strand=strand,mismatches=nm))
   status['alignments_checked']+=1
   if status['alignments_checked']%50000==0:save()
  assert proc.wait()==0,'Bowtie alignment failed'
 assert all(r['self_hits']==1 for r in queries.values()),'Expected exactly one self alignment per eligible query'
 for r in queries.values():r['no_other_locus_within_3_mismatches']=sum(r['off_target_counts'].values())==0
 (O/'results.json').write_text(json.dumps(dict(results=list(queries.values()),limitations=['Both strands of full ROS assembly, ungapped end-to-end Hamming distance <=3 only.','No PAM, nuclease, bulge, donor genotype, population alternate-genome off-target scan or editing safety approval.','PASS REF-span masks are conservative partial masks, not all indel effects.']),indent=2))
 status['summary']=[dict(window=w,queries=sum(r['window']==w for r in queries.values()),no_other_locus_within_3_mismatches=sum(r['window']==w and r['no_other_locus_within_3_mismatches'] for r in queries.values())) for w in sorted(masks)];status['state']='completed'
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:save();print(json.dumps(status))