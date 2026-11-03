from pathlib import Path
import json,gzip,urllib.request,hashlib,subprocess,math,traceback,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent;O=N/'ANO2_conservation_identity_20261009';O.mkdir(exist_ok=False);chain=Path('/mnt/d/Jin2024_work/zoonomia_hal/liftover/GCF_014441545.1ToCanFam3.over.chain.gz');start=39708074;end=39709074;mapping=[set() for _ in range(1000)]
with gzip.open(chain,'rt') as f:
 for line in f:
  x=line.split()
  if not x:continue
  if x[0]=='chain':h=x;tp=int(x[5]);qp=int(x[10]);assert x[4]=='+';continue
  length=int(x[0])
  if h[2]=='NC_051831.1':
   for pos in range(max(tp,start),min(tp+length,end)):
    q=qp+pos-tp;q=q if h[9]=='+' else int(h[8])-1-q;mapping[pos-start].add((h[7],q,h[9]))
  tp+=length;qp+=length
  if len(x)==3:tp+=int(x[1]);qp+=int(x[2])
assert all(len(v)==1 for v in mapping);m=[next(iter(v)) for v in mapping];assert all(v[0]==m[0][0] and v[2]==m[0][2] and v[1]==m[0][1]+i*(1 if m[0][2]=='+' else -1) for i,v in enumerate(m));chrom=m[0][0];lo=min(v[1] for v in m);hi=max(v[1] for v in m)+1
url=f'https://api.genome.ucsc.edu/getData/sequence?genome=canFam3;chrom={chrom};start={lo};end={hi}'
with urllib.request.urlopen(url,timeout=90) as r:raw=r.read()
(O/'canFam3.response.json').write_bytes(raw);dna=json.loads(raw)['dna'].upper();assert len(dna)==1000
source=Path('/mnt/d/Jin2024_work/zoonomia_hal/maf_dog/cand_36.chr27.39320764.maf');seen={};parts=[];buf=[]
def consume(lines):
 rows=[l.split() for l in lines if l.startswith('s ')];refs=[r for r in rows if r[1]=='Canis_lupus_familiaris.'+chrom]
 if not refs:return
 assert len(refs)==1;ref=refs[0];assert ref==rows[0] and ref[4]=='+';s=int(ref[2]);size=int(ref[3]);left=max(s,lo);right=min(s+size,hi)
 if left>=right:return
 assert len({len(r[6]) for r in rows})==1
 columns=[];pos=s
 for i,b in enumerate(ref[6]):
  if b!='-':
   if left<=pos<right:columns.append(i)
   pos+=1
 assert pos==s+size and columns
 a=columns[0];b=columns[-1]+1;out=[]
 for row in rows:
  seq=row[6][a:b];n=len(seq.replace('-',''))
  if n==0:continue
  offset=int(row[2])+len(row[6][:a].replace('-',''));out.append('s '+ ' '.join([row[1],str(offset),str(n),row[4],row[5],seq]))
 refdna=ref[6][a:b].replace('-','').upper();assert len(refdna)==right-left
 for i,bp in enumerate(refdna):assert left+i not in seen;seen[left+i]=bp
 parts.append('a\n'+'\n'.join(out)+'\n\n')
for line in source.open():
 if not line.strip():consume(buf);buf=[]
 elif not line.startswith('#'):buf.append(line.rstrip('\n'))
consume(buf);assert set(seen)==set(range(lo,hi));assert ''.join(seen[i] for i in range(lo,hi))==dna
maf=O/'ANO2_39708.maf';maf.write_text('##maf version=1\n\n'+''.join(parts));model=P/'06_v1.21.9/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod';wig=O/'ANO2_39708.wig';cmd=['/home/stardrako/miniforge3/envs/cactus/bin/phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(maf)];(O/'command.json').write_text(json.dumps(cmd))
with wig.open('w') as out,(O/'phyloP.err').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
vals={};pos=None
for line in wig.read_text().splitlines():
 if line.startswith('fixedStep'):f=dict(x.split('=') for x in line.split()[1:]);pos=int(f['start'])-1;assert int(f.get('step',1))==int(f.get('span',1))==1
 elif line.strip() and not line.startswith(('#','track')):v=float(line);assert pos is not None and pos not in vals and math.isfinite(v);vals[pos]=v;pos+=1
assert set(vals)==set(seen);values=[vals[i] for i in range(lo,hi)];max50=max(sum(values[i:i+50])/50 for i in range(951));result=dict(state='completed',ROS_window=['NC_051831.1',start,end],CanFam3_window=[chrom,lo,hi],strand=m[0][2],reference_coverage=1000,reference_matches_independent_UCSC=True,scored_bp=1000,max50bp_mean=max50,exceeds_legacy_6_5=max50>6.5,hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [chain,source,maf,wig,model]},independent_DNA_url=url,limitations=['Re-analysis of retained alignment and published model, not independent evolutionary data or biological validation.','Coordinate chain uniqueness checked here; prior reciprocal review assessed separately.']);(O/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))