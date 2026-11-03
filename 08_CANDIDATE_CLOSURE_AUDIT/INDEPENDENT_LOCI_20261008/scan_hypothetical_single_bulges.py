from pathlib import Path
import json,subprocess,gzip,pysam,hashlib,collections,datetime,traceback,fcntl
N=Path(__file__).resolve().parent;O=N/'hypothetical_NGG_single_bulge_20261009';O.mkdir(exist_ok=False)
lock=open('/tmp/safeharbor_NGG_bulges_20261009.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
source=N/'hypothetical_NGG_panel_20261009/results.json';rows=json.loads(source.read_text())['results'];guides={}
for i,g in enumerate(rows):
 if sum(g['other_NGG_counts'][str(k)]+g['unknown_PAM_counts'][str(k)] for k in range(3))==0:guides['g'+str(i)]=g
R=N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';fa=pysam.FastaFile(str(R));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
queries={};lookup={}
for gid,g in guides.items():
 seq=g['sequence'];templates=set(seq[:i]+seq[i+1:] for i in range(20))|set(seq[:i]+b+seq[i:] for i in range(21) for b in 'ACGT')
 for seq2 in sorted(templates):
  key='q'+str(len(queries));queries[key]=seq2;lookup[key]=gid
(O/'queries.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in queries.items()));(O/'guides.json').write_text(json.dumps(guides,indent=2));(O/'query_provenance.json').write_text(json.dumps(lookup));status=dict(state='running',hypothesis='Canonical NGG20nt panel; one DNA or RNA bulge plus <=2 substitutions',guides=len(guides),queries=len(queries),raw_alignments=0,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
def save():(O/'status.json').write_text(json.dumps(status,indent=2))
save();cmd=['/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin/bowtie','-f','-v','2','-a','--threads','4','/mnt/d/Jin2024_work/bowtie1_ROS_20261009/ROS',str(O/'queries.fa')];(O/'command.json').write_text(json.dumps(cmd));hits={gid:{} for gid in guides}
try:
 with (O/'scan.err').open('w') as err,gzip.open(O/'raw_alignments.tsv.gz','wt',compresslevel=1) as raw:
  proc=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True)
  for line in proc.stdout:
   raw.write(line);x=line.rstrip('\n').split('\t');gid=lookup[x[0]];g=guides[gid];s=int(x[3]);length=len(queries[x[0]]);strand=x[1];chrom=x[2];status['raw_alignments']+=1
   target=fa.fetch(chrom,s,s+length).upper();target=target if strand=='+' else rc(target)
   assert len(target)==length;assert sum(a!=b for a,b in zip(queries[x[0]],target))<=2
   if length==19:dist=min(sum(a!=b for a,b in zip(g['sequence'][:i]+g['sequence'][i+1:],target)) for i in range(20));kind='RNA_bulge'
   else:dist=min(sum(a!=b for a,b in zip(g['sequence'],target[:i]+target[i+1:])) for i in range(21));kind='DNA_bulge'
   assert dist<=2
   if chrom==g['chrom'] and s<g['start0']+20 and s+length>g['start0']:continue
   size=fa.get_reference_length(chrom)
   pp=fa.fetch(chrom,s+length,s+length+3).upper() if strand=='+' and s+length+3<=size else (rc(fa.fetch(chrom,s-3,s).upper()) if strand=='-' and s>=3 else None)
   if pp is not None and (pp[1] not in 'GN' or pp[2] not in 'GN'):continue
   certainty='NGG' if pp is not None and set(pp)<=set('ACGT') else 'unknown_PAM'
   hk=(chrom,s,length,strand);old=hits[gid].get(hk)
   if old is None or dist<old['substitutions']:hits[gid][hk]=dict(chrom=chrom,start0=s,end0=s+length,strand=strand,PAM=pp,PAM_status=certainty,bulge=kind,substitutions=dist)
   if status['raw_alignments']%50000==0:save()
  assert proc.wait()==0
 result=[]
 for gid,g in guides.items():
  hh=list(hits[gid].values());result.append(dict(guide_id=gid,window=g['window'],start0=g['start0'],strand=g['strand'],one_bulge_NGG_or_unknown_hits=hh,counts={str(k):sum(h['substitutions']==k for h in hh) for k in range(3)}))
 (O/'results.json').write_text(json.dumps(dict(results=result,limitations=['Hypothetical canonical NGG only, not selected nuclease or approved guides.','Exactly one base insertion/deletion and at most2 substitutions; no multi-base bulges or alternate genomes/PAMs.','Overlapping intended20nt locus omitted as intended-site alternative alignments; not evidence of absent local variation.','Single-gap distance independently recomputed from actual reference sequence; templates may report same alignment repeatedly, deduplicated.']),indent=2));status.update(state='completed',summary=[dict(window=w,guides=sum(r['window']==w for r in result),no_single_bulge_hit=sum(r['window']==w and not r['one_bulge_NGG_or_unknown_hits'] for r in result)) for w in sorted({r['window'] for r in result})])
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:save();print(json.dumps(status))
