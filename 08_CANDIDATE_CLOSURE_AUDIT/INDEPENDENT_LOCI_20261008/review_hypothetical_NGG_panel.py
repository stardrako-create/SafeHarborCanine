from pathlib import Path
import json,gzip,mmap,hashlib,collections,random,pysam,datetime
N=Path(__file__).resolve().parent; S=N/'mismatch_specificity_20261009'; O=N/'hypothetical_NGG_panel_20261009';O.mkdir(exist_ok=False)
assert json.loads((S/'status.json').read_text())['state']=='completed'
R=N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'
f=R.open('rb'); mm=mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ); idx={}
for line in Path(str(R)+'.fai').read_text().splitlines():
 x=line.split('\t');idx[x[0]]=tuple(map(int,x[1:5]))
def fetch(c,s,e):
 size,base,bases,width=idx[c]
 if s<0 or e>size:return None
 a=base+s//bases*width+s%bases;b=base+(e-1)//bases*width+(e-1)%bases+1
 seq=mm[a:b].replace(b'\n',b'').replace(b'\r',b'').decode().upper();assert len(seq)==e-s;return seq
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def pam(c,s,strand):
 seq=fetch(c,s+20,s+23) if strand=='+' else fetch(c,s-3,s)
 return seq if seq is None or strand=='+' else rc(seq)
def cls(seq):
 if seq is None:return 'unknown'
 if len(seq)!=3:return 'unknown'
 if any(b not in 'ACGTN' for b in seq):return 'unknown'
 if seq[1] not in 'GN' or seq[2] not in 'GN':return 'absent'
 return 'NGG' if 'N' not in seq else 'unknown'
rows=json.loads((S/'results.json').read_text())['results']; q={r['window']+'__'+str(r['offset0']):r for r in rows}
masks={r['window']:set(r['mask_offsets']) for r in json.loads((N/'two_loci_REF_mask_specificity.json').read_text())['results']}
a=N.parent/'ALL_CANDIDATES_20261005';ano='NC_051831.1_39708074_39709074';masks[ano]=set(next(r for r in json.loads((a/'masked_specificity_20261008.json').read_text())['results'] if r['id']==ano)['masked_offsets'])
fa=pysam.FastaFile(str(R));rng=random.Random(91)
for _ in range(200):
 c=rng.choice(list(idx));size=idx[c][0]
 if size<23:continue
 s=rng.randrange(size-22);assert fetch(c,s,s+23)==fa.fetch(c,s,s+23).upper()
guides={}
for k,r in q.items():
 for strand in '+-':
  po=r['offset0']+20 if strand=='+' else r['offset0']-3
  if po<0 or po+3>1000 or any(i in masks[r['window']] for i in range(po,po+3)):continue
  pp=pam(r['chrom'],r['start0'],strand)
  if cls(pp)!='NGG':continue
  guides[k+strand]=dict(query=k,window=r['window'],chrom=r['chrom'],start0=r['start0'],strand=strand,sequence=r['sequence'] if strand=='+' else rc(r['sequence']),PAM=pp,other_NGG_counts={str(i):0 for i in range(4)},unknown_PAM_counts={str(i):0 for i in range(4)},first10_other_NGG=[])
checked=0
with gzip.open(S/'alignments.tsv.gz','rt') as stream:
 for line in stream:
  x=line.rstrip('\n').split('\t');r=q[x[0]];c=x[2];s=int(x[3]);nm=0 if len(x)<8 or not x[7] else len(x[7].split(','))
  assert 0<=nm<=3
  if c==r['chrom'] and s==r['start0']:continue
  for orient in '+-':
   g=guides.get(x[0]+orient)
   if g is None:continue
   strand=x[1] if orient=='+' else ('-' if x[1]=='+' else '+'); pp=pam(c,s,strand);kind=cls(pp)
   if kind=='NGG':
    g['other_NGG_counts'][str(nm)]+=1
    if len(g['first10_other_NGG'])<10:g['first10_other_NGG'].append(dict(chrom=c,start0=s,strand=strand,PAM=pp,mismatches=nm))
   elif kind=='unknown':g['unknown_PAM_counts'][str(nm)]+=1
  checked+=1
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),hypothesis='20nt protospacer with canonical 3prime NGG PAM, both guide orientations; hypothetical design panel, not selected nuclease or approved guides',raw_alignment_sha256=hashlib.sha256((S/'alignments.tsv.gz').read_bytes()).hexdigest(),reference_sha256=json.loads((S/'status.json').read_text())['reference_sha256'],mmap_fetch_selftest=200,off_locus_alignments_checked=checked,results=list(guides.values()),summary=[],limits=['Ungapped Hamming <=3 on ROS reference only; no bulges or alternate PAM recognition.','No population alternate-genome off-target search or nuclease-specific activity score.','Unknown PAM sequence retained as uncertainty, never counted absent.','PAM and protospacer inside audited window and outside existing REF masks; not all population variants or haplotypes cleared.','No guide or biological safety approval; actual nuclease/cassette and empirical assays still needed.'])
for w in sorted(masks):
 rr=[g for g in guides.values() if g['window']==w];result['summary'].append(dict(window=w,hypothetical_NGG_guides=len(rr),no_other_NGG_or_unknown_PAM_within_mismatches={str(k):sum(sum(g['other_NGG_counts'][str(i)]+g['unknown_PAM_counts'][str(i)] for i in range(k+1))==0 for g in rr) for k in range(4)}))
(O/'results.json').write_text(json.dumps(result,indent=2));(O/'script.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest());print(json.dumps(result['summary']))
