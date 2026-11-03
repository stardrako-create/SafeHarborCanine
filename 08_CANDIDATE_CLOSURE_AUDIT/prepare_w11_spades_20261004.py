from pathlib import Path
import json,hashlib,datetime
A=Path(__file__).resolve().parent;O=A/'W11_SPADES_20261004';I=A/'SV_PBGV000012/ROS_expanded_junction_20261001/reads.fastq';pairs={}
with I.open() as f:
 while True:
  h=f.readline().strip()
  if not h:break
  seq=f.readline().strip();plus=f.readline().strip();q=f.readline().strip();assert plus=='+' and len(seq)==len(q)
  name,m=h[1:].rsplit('_mate',1);assert m in ['1','2'] and m not in pairs.get(name,{})
  pairs.setdefault(name,{})[m]=(seq,q)
n=single=0
with (O/'R1.fastq').open('x') as a,(O/'R2.fastq').open('x') as b,(O/'single.fastq').open('x') as c:
 for name,reads in sorted(pairs.items()):
  if len(reads)==2:
   n+=1
   for m,f in [('1',a),('2',b)]:
    seq,q=reads[m];f.write('@'+name+'/'+m+'\n'+seq+'\n+\n'+q+'\n')
  else:
   single+=1;m=next(iter(reads));seq,q=reads[m];c.write('@'+name+'/'+m+'\n'+seq+'\n+\n'+q+'\n')
assert 2*n+single==14716
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sample':'PBGV000012','pairs':n,'single_reads':single,'total_reads':2*n+single,'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [I,O/'R1.fastq',O/'R2.fastq',O/'single.fastq',Path(__file__)]},'status':'inputs_prepared_assembler_download_timeout','limits':['Previously selected local windows and retained mates only; not full raw read recruitment']};(O/'input_receipt.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))