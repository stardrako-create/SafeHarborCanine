from pathlib import Path
import json,random,subprocess,hashlib
N=Path(__file__).resolve().parent;O=N/'hypothetical_TTTV23_panel_20261009/toy_validation';O.mkdir(exist_ok=False);rng=random.Random(419);ref=''.join(rng.choice('ACGT') for _ in range(160));ref+=ref[20:50];rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1];B='/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin/';(O/'ref.fa').write_text('>toy\n'+ref+'\n');cmd=[B+'bowtie-build','--quiet',str(O/'ref.fa'),str(O/'index')];subprocess.run(cmd,check=True,capture_output=True,timeout=120);results=[]
for lengths,v in [([22,24],2),([23],3)]:
 q={}
 for length in lengths:
  for start in [0,20,80]:
   for strand,seq in [('+',ref[start:start+length]),('-',rc(ref[start:start+length]))]:
    for nm in range(v+1):
     t=list(seq)
     for z in range(nm):t[z]='C' if t[z]!='C' else 'A'
     q[f'{length}_{start}_{strand}_{nm}']=''.join(t)
 stem='lengths_'+str(lengths[0]);f=O/(stem+'.fa');f.write_text(''.join('>'+k+'\n'+s+'\n' for k,s in q.items()));cmd=[B+'bowtie','-f','-v',str(v),'-a','--quiet',str(O/'index'),str(f)];ret=subprocess.run(cmd,check=True,capture_output=True,text=True,timeout=120);(O/(stem+'.command.json')).write_text(json.dumps(cmd));(O/(stem+'.out')).write_text(ret.stdout);(O/(stem+'.err')).write_text(ret.stderr)
 actual={(x[0],x[1],int(x[3])) for line in ret.stdout.splitlines() for x in [line.split('\t')]};expected={(k,strand,s) for k,seq in q.items() for strand,ori in [('+',seq),('-',rc(seq))] for s in range(len(ref)-len(seq)+1) if sum(a!=b for a,b in zip(ori,ref[s:s+len(seq)]))<=v};assert actual==expected;results.append(dict(lengths=lengths,mismatch_limit=v,queries=len(q),alignments=len(actual),matches_bruteforce=True))
(N/'hypothetical_TTTV23_search_validation_20261009.json').write_text(json.dumps(dict(state='completed',results=results,limits=['Controlled search-completeness test, not empirical nuclease specificity or real genome uniqueness validation.']),indent=2));print(json.dumps(results))
