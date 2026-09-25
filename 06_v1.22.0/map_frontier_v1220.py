from pathlib import Path
import gzip,json,csv,hashlib
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.0';O.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
front=read(P/'06_v1.21.9/local_tradeoff_frontier_annotated_v1219.tsv');regs={r['left_gene']:r for r in read(P/'06_v1.21.8/candidates_scored_v1218.tsv') if r['left_gene'] in ['ANO2','LOC119876429','NPNT']}
chains=json.loads((P/'06_v1.21.4/chain_blocks.json').read_text());first={g:max(chains[g],key=lambda c:c['covered']) for g in regs}
def chainblocks(path):
 with gzip.open(path,'rt') as f:
  for line in f:
   v=line.split()
   if not v:continue
   if v[0]=='chain':h=v;tp=int(v[5]);qp=int(v[10]);continue
   n=int(v[0]);assert h[4]=='+'
   qa,qb=(qp,qp+n) if h[9]=='+' else (int(h[8])-qp-n,int(h[8])-qp)
   yield h,tp,tp+n,qa,qb
   if len(v)==3:tp+=n+int(v[1]);qp+=n+int(v[2])
# Two independently composed coordinate routes; retain all mappings to reject ambiguity.
ros3={};qneeded={}
for g,r in regs.items():
 c=first[g];assert c['strand']=='+';d={}
 for a,b,x,y in c['blocks']:
  for t,q in zip(range(a,b),range(x,y)):
   if int(r['start'])<=t<int(r['end']):d[t]=q
 ros3[g]=d;qneeded[g]=set(d.values())
forward={g:{} for g in regs};direct={g:{} for g in regs}
bounds={g:(min(v),max(v)+1) for g,v in qneeded.items()}
cf3chrom={'ANO2':'chr27','LOC119876429':'chr7','NPNT':'chr32'}
for h,a,b,x,y in chainblocks(P/'06_v1.21.7/canFam3ToCanFam6.over.chain.gz'):
 for g in regs:
  if h[2]!=cf3chrom[g]:continue
  needed=qneeded[g]
  low,high=bounds[g]
  if b<=low or a>=high:continue
  for q in range(max(a,low),min(b,high)):
   if q in needed:forward[g].setdefault(q,set()).add((h[7],x+q-a if h[9]=='+' else y-1-(q-a)))
for h,a,b,x,y in chainblocks(P/'06_v1.21.7/canFam6ToGCF_014441545.1.over.chain.gz'):
 for g,r in regs.items():
  if h[7]!=r['chrom']:continue
  for t in range(max(x,int(r['start'])),min(y,int(r['end']))):direct[g].setdefault(t,set()).add((h[2],a+t-x if h[9]=='+' else a+y-1-t))
out=[];targets=[];partial_maps=[]
for r in front:
 g=r['genes'].split('/')[0];s=int(r['start']);e=int(r['end']);dest=[];ok=True
 for t in range(s,e):
  f=forward[g].get(ros3[g].get(t),set());d=direct[g].get(t,set())
  if len(f)!=1 or f!=d:ok=False;continue
  dest.append(next(iter(f)))
 contiguous=ok and len({c for c,p in dest})==1 and len({p for c,p in dest})==e-s and max(p for c,p in dest)-min(p for c,p in dest)+1==e-s
 r.update(canfam6_mapping_status='two_routes_unique_identical_contiguous' if contiguous else 'incomplete_ambiguous_or_noncontiguous_not_quantified',canfam6_chrom=dest[0][0] if contiguous else '',canfam6_start=min(p for c,p in dest) if contiguous else '',canfam6_end=max(p for c,p in dest)+1 if contiguous else '')
 blocks=[]
 for chrom,pos in sorted(set(dest)):
  if blocks and blocks[-1][0]==chrom and blocks[-1][2]==pos:blocks[-1][2]=pos+1
  else:blocks.append([chrom,pos,pos+1])
 partial_maps.append(dict(genes=r['genes'],chrom=r['chrom'],start=s,end=e,shared_unique_destination_bp=len(set(dest)),source_bp_with_two_route_agreement=len(dest),source_bp=e-s,blocks=blocks))
 out.append(r)
aliases={v[0]:v[3] for line in (P/'06_v1.21.7/canFam6.chromAlias.txt').read_text().splitlines() if not line.startswith('#') and len(v:=line.split('\t'))>=4}
for g in regs:
 bb=[b for r in partial_maps if r['genes'].split('/')[0]==g for b in r['blocks']];assert len({b[0] for b in bb})==1
 if bb:targets.append(dict(name=g+'_frontier_span',kind='exploratory_frontier_span',chrom=aliases[bb[0][0]],start=min(b[1] for b in bb),end=max(b[2] for b in bb)))
with (O/'frontier_canfam6_mapping.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
(O/'fragment_targets.json').write_text(json.dumps(targets,indent=2));summary=dict(windows=len(out),accepted=sum(r['canfam6_mapping_status']=='two_routes_unique_identical_contiguous' for r in out),targets=targets,chain_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (P/'06_v1.21.7').glob('*.chain.gz')},limitations=['Exact agreement between two chain routes; no sequence or individual donor genotype confirmation','Unmapped/ambiguous/noncontiguous windows remain unquantified, not zero'])
(O/'frontier_mapping_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
(O/'frontier_shared_mapping_blocks.json').write_text(json.dumps(partial_maps,indent=2))
