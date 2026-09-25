from pathlib import Path
import gzip,csv,json,hashlib
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.21.7'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
windows=read(P/'06_v1.21.5/local_canfam3_blocks.tsv');ros={r['genes']:r for r in read(P/'06_v1.21.5/cohort_and_support.tsv') if not r['genes'].endswith('_region')}
def project(path,inverse=False):
 results={w['genes']:{} for w in windows};h=None
 with gzip.open(path,'rt') as f:
  for line in f:
   v=line.split()
   if not v:continue
   if v[0]=='chain':
    h=v;tp=int(v[5]);qp=int(v[10]);continue
   size=int(v[0]);assert h[4]=='+'
   qa,qb=(qp,qp+size) if h[9]=='+' else (int(h[8])-qp-size,int(h[8])-qp)
   for w in windows:
    r=ros[w['genes']] if inverse else w;c=h[7] if inverse else h[2];a,b=(qa,qb) if inverse else (tp,tp+size)
    if c!=r['chrom']:continue
    s=max(a,int(r['start']));e=min(b,int(r['end']))
    if s>=e:continue
    if h[9]=='+':x,y=((tp+s-qa,tp+e-qa) if inverse else (qa+s-tp,qa+e-tp))
    else:x,y=((tp+qb-e,tp+qb-s) if inverse else (qb-(e-tp),qb-(s-tp)))
    d=results[w['genes']].setdefault(h[12],dict(chain_id=h[12],score=int(h[1]),chrom=h[2] if inverse else h[7],strand=h[9],mapped_bp=0,blocks=[]));d['mapped_bp']+=e-s;d['blocks'].append([s,e,x,y])
   if len(v)==3:tp+=size+int(v[1]);qp+=size+int(v[2])
 return results
forward=project(O/'canFam3ToCanFam6.over.chain.gz');direct=project(O/'canFam6ToGCF_014441545.1.over.chain.gz',True);summary=[]
for w in windows:
 g=w['genes'];f=max(forward[g].values(),key=lambda x:x['mapped_bp']);d=max(direct[g].values(),key=lambda x:x['mapped_bp'])
 fs=set(i for x in f['blocks'] for i in range(x[2],x[3]));ds=set(i for x in d['blocks'] for i in range(x[2],x[3]));same=f['chrom']==d['chrom'] and fs==ds
 assert same and f['mapped_bp']==d['mapped_bp']==1000,(g,f,d)
 assert max(fs)-min(fs)+1==1000
 summary.append(dict(genes=g,chrom=f['chrom'],start=min(fs),end=max(fs)+1,canfam3_chain_id=f['chain_id'],direct_chain_id=d['chain_id'],two_paths_identical=same,mapped_bp=len(fs)))
with (O/'local_canfam6_windows.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');wr.writeheader();wr.writerows(summary)
(O/'mapping_provenance.json').write_text(json.dumps(dict(forward=forward,direct_inverse=direct,summary=summary,chain_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.glob('*.chain.gz')},limitations='Coordinate agreement from two chain paths. No individual donor genotype or sequence-level confirmation in this donor.'),indent=2))
print(json.dumps(summary,indent=2))
