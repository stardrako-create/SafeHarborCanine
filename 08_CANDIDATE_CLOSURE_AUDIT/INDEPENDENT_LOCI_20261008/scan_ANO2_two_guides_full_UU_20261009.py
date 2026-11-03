from pathlib import Path
import subprocess,json,hashlib,datetime,fcntl,gzip,pysam,traceback,os
N=Path(__file__).resolve().parent;O=N/'ANO2_two_guides_full_UU_20261009';O.mkdir(exist_ok=False)
lock=open('/tmp/safeharbor_full_UU_two_guides.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
U=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');B=Path('/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin');I=Path('/mnt/d/Jin2024_work/bowtie1_UU_20261009');I.mkdir(exist_ok=False);idx=I/'UU'
status=dict(state='running',pid=os.getpid(),stage='building_UU_index',utc_start=datetime.datetime.now(datetime.timezone.utc).isoformat())
def save():
 p=O/'status.tmp';p.write_text(json.dumps(status,indent=2));p.replace(O/'status.json')
def run(cmd,name):
 (O/(name+'.command.json')).write_text(json.dumps(cmd))
 with (O/(name+'.out')).open('w') as out,(O/(name+'.err')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True)
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
try:
 save();run([str(B/'bowtie-build'),'--threads','4',str(U),str(idx)],'build');fa=pysam.FastaFile(str(U));guides=json.loads((N/'ANO2_TTTV23_REFERENCE_MODEL_CANDIDATES_20261009.json').read_text())['results'];hits={g['hypothesis_id']:{'ungapped':{},'bulge':{}} for g in guides};selfhits={g['hypothesis_id']:0 for g in guides}
 for g in guides:
  s=6741835+1000-(g['protospacer_start0']-39708074+23);g['UU_expected']=dict(chrom='NC_049248.1',start0=s,end0=s+23,strand='-')
  assert rc(fa.fetch('NC_049248.1',s,s+23).upper())==g['protospacer_5to3']
  assert rc(fa.fetch('NC_049248.1',s+23,s+27).upper())==g['PAM_5to3']
 def pamclass(p):
  if p is None:return 'unknown'
  if any(all(x=='N' or x in ('ACG' if y=='V' else y) for x,y in zip(p,pat)) for pat in ('TTTV','CTTV','TCTV','TTCV')):return 'unknown' if 'N' in p else 'compatible'
  return 'absent'
 for mode in ('ungapped','bulge'):
  status['stage']='search_'+mode;save();q={};look={}
  for g in guides:
   seq=g['protospacer_5to3'];templates={seq} if mode=='ungapped' else {seq[:i]+seq[i+1:] for i in range(23)}|{seq[:i]+b+seq[i:] for i in range(24) for b in 'ACGT'}
   for t in sorted(templates):key='q'+str(len(q));q[key]=t;look[key]=g
  qf=O/(mode+'.fa');qf.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in q.items()));cmd=[str(B/'bowtie'),'-f','-v','3' if mode=='ungapped' else '2','-a','--threads','4',str(idx),str(qf)];(O/(mode+'.command.json')).write_text(json.dumps(cmd));rawcount=0
  with gzip.open(O/(mode+'.tsv.gz'),'wt') as raw,(O/(mode+'.err')).open('w') as err:
   p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True)
   for line in p.stdout:
    raw.write(line);rawcount+=1;x=line.rstrip().split('\t');g=look[x[0]];e=g['UU_expected'];c=x[2];s=int(x[3]);strand=x[1];length=len(q[x[0]]);t=fa.fetch(c,s,s+length).upper();t=t if strand=='+' else rc(t);assert sum(a!=b for a,b in zip(q[x[0]],t))<=(3 if mode=='ungapped' else 2)
    if mode=='ungapped' and (c,s,strand)==(e['chrom'],e['start0'],e['strand']):assert t==g['protospacer_5to3'];selfhits[g['hypothesis_id']]+=1;continue
    if mode=='bulge' and c==e['chrom'] and s<e['end0'] and s+length>e['start0']:continue
    seq=g['protospacer_5to3'];nm=sum(a!=b for a,b in zip(seq,t)) if mode=='ungapped' else (min(sum(a!=b for a,b in zip(seq[:i]+seq[i+1:],t)) for i in range(23)) if length==22 else min(sum(a!=b for a,b in zip(seq,t[:i]+t[i+1:])) for i in range(24)))
    size=fa.get_reference_length(c);pp=fa.fetch(c,s-4,s).upper() if strand=='+' and s>=4 else (rc(fa.fetch(c,s+length,s+length+4).upper()) if strand=='-' and s+length+4<=size else None);cls=pamclass(pp)
    if cls=='absent':continue
    key=f'{c}:{s}:{strand}:{length}';d=hits[g['hypothesis_id']][mode];old=d.get(key)
    if old is None or nm<old['mismatches']:d[key]=dict(chrom=c,start0=s,end0=s+length,strand=strand,PAM=pp,PAM_class=cls,mismatches=nm)
   assert p.wait()==0
  status[mode+'_raw_alignments']=rawcount;save()
 assert all(n==1 for n in selfhits.values())
 results=[dict(guide=g,exact_intended_hits=selfhits[g['hypothesis_id']],ungapped_hits=list(hits[g['hypothesis_id']]['ungapped'].values()),single_bulge_hits=list(hits[g['hypothesis_id']]['bulge'].values())) for g in guides]
 (O/'results.json').write_text(json.dumps(dict(results=results,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),index_files=[dict(name=p.name,size=p.stat().st_size) for p in I.iterdir()],limits=['Complete UU reference search under defined exploratory PAM/mismatch/bulge model; actual nuclease not selected.','Intended interval overlapping bulge alignments excluded, not shifted local cleavage clearance.','Unknown PAM retained conservatively; N/gaps/missing assembly and population alleles not resolved.','Another reference assembly is not experimental independent validation.']),indent=2));status.update(state='completed',stage='review_complete',summary=[dict(guide=r['guide']['hypothesis_id'],ungapped_hits=len(r['ungapped_hits']),single_bulge_hits=len(r['single_bulge_hits'])) for r in results])
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:status['utc_end']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
