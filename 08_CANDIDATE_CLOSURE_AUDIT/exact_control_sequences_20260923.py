from pathlib import Path
import sys,json,csv,collections,time,fcntl,hashlib
sys.path.insert(0,'/home/stardrako/safeharbor_audit_dependencies')
import ahocorasick,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'08_CANDIDATE_CLOSURE_AUDIT/CONTROL_LOCUS_REVIEW_20260922/exact_sequences';O.mkdir(parents=True,exist_ok=True)
lock=open('/tmp/safeharbor_controls_exact_20260923.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
R=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';fa=pysam.FastaFile(str(R));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
rows=list(csv.DictReader((O.parent/'transcript_and_boundary_review.tsv').open(),delimiter='\t'));rows=[dict(r,window_id=r['id'],start=r['start0'],end=r['end0']) for r in rows];patterns={};windows=[]
for r in rows:
 s=int(r['start']);seq=fa.fetch(r['chrom'],s,int(r['end'])).upper()
 for k in [20,25,50,100]:
  for offset in range(1000-k+1):
   q=seq[offset:offset+k];assert set(q)<=set('ACGT');canonical=min(q,rc(q));patterns.setdefault(canonical,len(patterns));windows.append(dict(window_id=r['window_id'],chrom=r['chrom'],start0=s+offset,end0=s+offset+k,length=k,pattern_id=patterns[canonical]))
A=ahocorasick.Automaton()
for seq,i in patterns.items():
 for q in {seq,rc(seq)}:A.add_word(q,(i,len(q)))
A.make_automaton();counts=collections.Counter();locations=collections.defaultdict(list);MAX=100;BLOCK=1000000
def scan(sequence,blocksize,automaton):
 for start in range(0,len(sequence),blocksize):
  lo=max(0,start-MAX+1);chunk=sequence[lo:start+blocksize]
  for end,(i,k) in automaton.iter(chunk):
   absolute_end=lo+end+1
   if absolute_end>start:yield i,absolute_end-k,absolute_end
# Boundary/overlapping-hit self-test against brute force on both strands.
T=ahocorasick.Automaton();toy={'AAAA':0,'ACGT':1,'TTTT':0}
for q,i in toy.items():T.add_word(q,(i,len(q)))
T.make_automaton();toyseq='AAAAAACGTTTTTACGT';observed=collections.Counter(scan(toyseq,5,T));expected=collections.Counter((i,pos,pos+len(q)) for q,i in toy.items() for pos in range(len(toyseq)-len(q)+1) if toyseq[pos:pos+len(q)]==q);assert observed==expected
started=time.time();bases=0;done=[]
def status(state):
 tmp=O/'exact_scan_status.tmp';tmp.write_text(json.dumps(dict(status=state,started=started,updated=time.time(),bases_scanned=bases,contigs_completed=done,patterns=len(patterns),tested_segments=len(windows)),indent=2));tmp.replace(O/'exact_scan_status.json')
status('running')
for chrom in fa.references:
 length=fa.get_reference_length(chrom)
 for start in range(0,length,BLOCK):
  lo=max(0,start-MAX+1);seq=fa.fetch(chrom,lo,min(length,start+BLOCK)).upper()
  for end,(i,k) in A.iter(seq):
   absolute_end=lo+end+1
   if absolute_end<=start:continue
   counts[i]+=1
   if len(locations[i])<10:locations[i].append([chrom,absolute_end-k,absolute_end])
  bases+=min(BLOCK,length-start)
 done.append(chrom);status('running')
for r in windows:
 r['exact_locus_count_both_strands']=counts[r['pattern_id']];assert r['exact_locus_count_both_strands']>=1;r['unique_exact']=r['exact_locus_count_both_strands']==1
with (O/'all_short_segments_exact_uniqueness.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(windows[0]),delimiter='\t');w.writeheader();w.writerows(windows)
summary=[]
for r in rows:
 for k in [20,25,50,100]:
  vv=[x for x in windows if x['window_id']==r['window_id'] and x['length']==k];summary.append(dict(window_id=r['window_id'],length=k,segments=len(vv),unique=sum(x['unique_exact'] for x in vv),nonunique=sum(not x['unique_exact'] for x in vv),max_locus_count=max(x['exact_locus_count_both_strands'] for x in vv)))
(O/'exact_nonunique_locations.json').write_text(json.dumps({str(i):dict(count=n,first10_locations=locations[i]) for i,n in counts.items() if n>1},indent=2))
(O/'exact_sequence_review.json').write_text(json.dumps(dict(status='completed',summary=summary,bases_scanned=bases,contigs=len(done),seconds=time.time()-started,self_test='chunk boundaries and overlapping both-strand counts agree with brute force',reference=str(R),limitations=['Exact sequence matches only, not mismatch/bulge off-target analysis','Includes assembly contigs; alternate/haplotype sequences may inflate locus counts','No PAM, nuclease, guide or donor genotype specified','1000bp target windows do not define an insertion base','Unique exact sequence does not establish editing safety']),indent=2));status('completed');print(json.dumps(summary,indent=2))

