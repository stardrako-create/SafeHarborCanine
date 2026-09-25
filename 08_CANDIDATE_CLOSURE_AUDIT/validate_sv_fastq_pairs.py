"""Full paired FASTQ structure/identity check; never runs on partial downloads."""
from pathlib import Path
import argparse,gzip,json,time,itertools,os,msvcrt

def records(path):
 with gzip.open(path,'rb') as f:
  while True:
   h=f.readline()
   if not h:return
   s=f.readline().rstrip(b'\r\n');plus=f.readline();q=f.readline().rstrip(b'\r\n')
   if not h.startswith(b'@') or not plus.startswith(b'+') or not s or len(s)!=len(q):
    raise ValueError('Invalid FASTQ record in '+str(path))
   if min(q)<33 or max(q)>126:raise ValueError('Invalid quality encoding')
   ident=h[1:].split()[0]
   if ident.endswith((b'/1',b'/2')):ident=ident[:-2]
   yield ident,s,q

def validate(a,b,progress=None):
 counts={'pairs':0,'R1_bases':0,'R2_bases':0,'R1_N':0,'R2_N':0,'R1_Q30':0,'R2_Q30':0}
 for x,y in itertools.zip_longest(records(a),records(b)):
  if x is None or y is None:raise ValueError('Different mate record counts')
  if x[0]!=y[0]:raise ValueError('Mate identifiers do not match at pair '+str(counts['pairs']+1))
  counts['pairs']+=1
  for prefix,r in [('R1',x),('R2',y)]:
   counts[prefix+'_bases']+=len(r[1]);counts[prefix+'_N']+=r[1].upper().count(b'N');counts[prefix+'_Q30']+=len(r[2].translate(None,bytes(range(63))))
  if progress and counts['pairs']%100000==0:progress(counts)
 if not counts['pairs']:raise ValueError('Empty FASTQ inputs')
 return counts

def main():
 ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);args=ap.parse_args();d=args.directory
 status=json.loads((d/'status.json').read_text());assert status['status']=='download_complete'
 lock=(d/'pair_validation.lock').open('a+b');lock.seek(0);msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
 def progress(counts,state='running'):
  tmp=d/'pair_validation_status.tmp';tmp.write_text(json.dumps(dict(status=state,pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**counts),indent=2));os.replace(tmp,d/'pair_validation_status.json')
 progress({'pairs':0})
 paths=[d/f'SRR15734832_{i}.fastq.gz' for i in [1,2]]
 for p in paths:
  receipt=json.loads((d/(p.name+'.validated.json')).read_text())
  assert receipt['validated'] and p.stat().st_size==receipt['bytes']
 try:out=validate(*paths,progress=progress)
 except Exception as e:
  progress({'error':repr(e)},'failed');raise
 progress(out,'completed')
 out.update(status='validated_FASTQ_structure_and_mate_identity',gzip_streams_read_to_EOF=True,completed_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),scope='Input QC only; no alignment, sample identity or SV validation')
 (d/'fastq_pair_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))

if __name__=='__main__':main()
