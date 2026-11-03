from pathlib import Path
import json,subprocess,random,hashlib
N=Path(__file__).resolve().parent;O=N/'hypothetical_NGG_single_bulge_20261009';S=O/'results.json';status=json.loads((O/'status.json').read_text());assert status['state']=='completed'
rng=random.Random(349);dna=''.join(rng.choice('ACGT') for _ in range(150));rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1];q={}
for l in [19,21]:
 for off in [0,25,60]:
  seq=dna[off:off+l]
  for strand,ss in [('+',seq),('-',rc(seq))]:
   for nm in [0,1,2]:
    t=list(ss)
    for i in range(nm):t[i]='A' if t[i]!='A' else 'C'
    q[f'{l}_{off}_{strand}_{nm}']=''.join(t)
toy=O/'toy_validation';toy.mkdir(exist_ok=False);(toy/'ref.fa').write_text('>toy\n'+dna+'\n');(toy/'q.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in q.items()));B='/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin/'
commands=[[B+'bowtie-build','--quiet',str(toy/'ref.fa'),str(toy/'index')],[B+'bowtie','-f','-v','2','-a','--quiet',str(toy/'index'),str(toy/'q.fa')]]
for i,cmd in enumerate(commands):
 (toy/f'{i}.command.json').write_text(json.dumps(cmd));ret=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=120);(toy/f'{i}.out').write_text(ret.stdout);(toy/f'{i}.err').write_text(ret.stderr)
actual={(x[0],x[1],int(x[3])) for line in (toy/'1.out').read_text().splitlines() for x in [line.split('\t')]}
expected={(k,strand,s) for k,seq in q.items() for strand,ori in [('+',seq),('-',rc(seq))] for s in range(len(dna)-len(seq)+1) if sum(a!=b for a,b in zip(ori,dna[s:s+len(seq)]))<=2};assert actual==expected
rr=json.loads(S.read_text())['results'];summaries=[]
for w in sorted({r['window'] for r in rr}):
 subset=[r for r in rr if r['window']==w];summaries.append(dict(window=w,guides=len(subset),no_other_locus_single_bulge_NGG_or_unknown_through_substitutions={str(k):sum(sum(r['counts'][str(i)] for i in range(k+1))==0 for r in subset) for k in range(3)},total_deduplicated_hits=sum(len(r['one_bulge_NGG_or_unknown_hits']) for r in subset)))
result=dict(state='completed',selftest=dict(query_lengths=[19,21],queries=len(q),alignments=len(actual),matches_exhaustive_bruteforce=True),summary=summaries,source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),limits=['Hypothetical canonicalNGG reference-only one-base-bulge analysis, not guide approval or actual nuclease choice.','Counts exclude alignments overlapping intended20nt interval, not a proof of absence of local shifted cleavage.','No new locus veto or zero-offtarget requirement imposed.']);out=N/'hypothetical_NGG_single_bulge_interpretation_20261009.json';assert not out.exists();out.write_text(json.dumps(result,indent=2));print(json.dumps(result))
