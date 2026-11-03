from pathlib import Path
import json,hashlib,pysam,datetime
N=Path(__file__).resolve().parent
S=N/'ANO2_TTTV23_offsite_UU_and_catalogue_20261009/results.json'
O=N/'ANO2_offsite_three_indels_20261009'; O.mkdir(exist_ok=False)
f=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa')
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
pams={p+b for p in ('TTT','CTT','TCT','TTC') for b in 'ACG'}
def distance(g,t):
 if len(t)==23:return sum(a!=b for a,b in zip(g,t))
 if len(t)==22:return min(sum(a!=b for a,b in zip(g[:i]+g[i+1:],t)) for i in range(23))
 return min(sum(a!=b for a,b in zip(g,t[:i]+t[i+1:])) for i in range(24))
def scan(seq,g):
 hits=[]
 for strand,s in [('+',seq),('-',rc(seq))]:
  for i in range(len(s)-25):
   if s[i:i+4] not in pams:continue
   for length in (22,23,24):
    t=s[i+4:i+4+length]
    if len(t)!=length or set(t)-set('ACGT'):continue
    d=distance(g,t)
    if d<=(3 if length==23 else 2):hits.append(dict(strand=strand,oriented_offset=i,PAM=s[i:i+4],target=t,length=length,mismatches=d))
 return hits
out=[]
for r in json.loads(S.read_text())['results']:
 for q in r.get('catalogue',[]):
  for a in q.get('alleles',[]):
   if a['status']=='single_SNP_model_assessed':continue
   c=r['UU_core']['chrom']; pos=a['pos1']-1; lo=max(0,pos-100); hi=pos+len(a['ref'])+100
   ref=f.fetch(c,lo,hi).upper(); off=pos-lo
   assert ref[off:off+len(a['ref'])]==a['ref'].upper()
   assert set(a['alt'].upper())<=set('ACGT')
   alt=ref[:off]+a['alt'].upper()+ref[off+len(a['ref']):]
   g=r['guide']['protospacer_5to3']; before=scan(ref,g); after=scan(alt,g)
   (O/(r['site_id']+'.fa')).write_text('>UU_reference\n'+ref+'\n>single_catalogue_ALT\n'+alt+'\n')
   out.append(dict(site=r['site_id'],guide=r['guide']['hypothesis_id'],chrom=c,context_start0=lo,context_end0=hi,allele=a,reference_hits=before,alternate_hits=after,ref_sha256=hashlib.sha256(ref.encode()).hexdigest(),alt_sha256=hashlib.sha256(alt.encode()).hexdigest()))
assert len(out)==3
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),results=out,limits=['Three queried catalogued indels applied individually to UU reference contexts, no phase or combined haplotype inference.','Both strands: defined PAM family and ungapped 23nt <=3 substitutions or one-base bulge <=2 substitutions.','This closes only these three allele contexts; no whole-population or all alternate-genome clearance and no empirical cleavage claim.'])
(O/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(state=result['state'],sites=[dict(site=x['site'],reference_hits=len(x['reference_hits']),alternate_hits=len(x['alternate_hits'])) for x in out])))
