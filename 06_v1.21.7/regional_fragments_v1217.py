from pathlib import Path
import urllib.request,gzip,struct,json,csv,hashlib,concurrent.futures,threading
O=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino\06_v1.21.7');C=O/'range_cache';C.mkdir(exist_ok=True)
URL='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7807nnn/GSM7807442/suppl/GSM7807442_atac_fragments.tsv.gz';SIZE=1920882177;MAGIC=b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00'
partial=O/'GSM7807442_atac_fragments.tsv.gz.partial';requests=[];lock=threading.Lock()
def fetch(a,b):
 b=min(b,SIZE-1);p=C/f'{a}-{b}.bgzfpart'
 if p.exists():return p.read_bytes()
 if partial.exists() and b<partial.stat().st_size:
  with partial.open('rb') as f:f.seek(a);d=f.read(b-a+1)
 else:
  with urllib.request.urlopen(urllib.request.Request(URL,headers={'Range':f'bytes={a}-{b}'}),timeout=45) as r:
   assert r.status==206 and r.headers['Content-Range']==f'bytes {a}-{b}/{SIZE}'
   d=r.read(b-a+1);assert len(d)==b-a+1
 p.write_bytes(d);return d
def decode(d):
 pos=d.find(MAGIC);assert pos>=0;parts=[]
 while pos+18<=len(d) and d[pos:pos+16]==MAGIC:
  n=struct.unpack_from('<H',d,pos+16)[0]+1
  if pos+n>len(d):break
  parts.append(gzip.decompress(d[pos:pos+n]));pos+=n
 lines=b''.join(parts).decode().split('\n')[1:-1]
 return [l.split('\t') for l in lines if l and not l.startswith('#')]
head=gzip.decompress(fetch(0,struct.unpack_from('<H',partial.open('rb').read(18),16)[0])).decode()
contigs=[l.split('=',1)[1] for l in head.splitlines() if l.startswith('# primary_contig=')];rank={c:i for i,c in enumerate(contigs)}
def key(r):return rank.get(r[0],len(rank)),int(r[1])
def sample(pos):
 rows=decode(fetch(pos,pos+131071));assert rows
 return key(rows[0])
def query(t):
 target=(rank[t['chrom']],max(0,t['start']-10000));lo=0;hi=SIZE-131072;steps=[]
 while hi-lo>32768:
  mid=(lo+hi)//2;k=sample(mid);steps.append([mid,list(k)])
  if k<target:lo=mid
  else:hi=mid
 a=max(0,lo-131072);b=min(SIZE-1,hi+1048576);rows=decode(fetch(a,b));keys=[key(r) for r in rows];assert keys==sorted(keys),'Local coordinate order mismatch'
 assert keys[0]<(rank[t['chrom']],t['start']) and keys[-1]>(rank[t['chrom']],t['end']),'Window not bracketed'
 selected=[r for r in rows if r[0]==t['chrom'] and int(r[1])<t['end'] and int(r[2])>t['start']]
 p=O/(t['name'].replace('/','_')+'.fragments.tsv');p.write_text('\n'.join('\t'.join(r) for r in selected)+'\n')
 result=dict(target=t,search_steps=steps,range_start=a,range_end=b,first_key=list(keys[0]),last_key=list(keys[-1]),fragment_rows=len(selected),starts_in_window=sum(t['start']<=int(r[1])<t['end'] for r in selected),raw_file=p.name,interpretation='Complete starts-in-window query under coordinate-sorted file assumption; overlap counts limited to fetched bracket, so extremely long fragments starting before it could be omitted.')
 (O/(p.stem+'.query.json')).write_text(json.dumps(result,indent=2));print(t['name'],len(selected),result['starts_in_window'],flush=True);return result
targets=json.loads((O/'fragment_targets.json').read_text())
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:result=list(pool.map(query,targets))
(O/'regional_fragment_queries.json').write_text(json.dumps(dict(url=URL,total_remote_bytes=SIZE,contig_order_from_header=contigs,queries=result,cache={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in C.iterdir()}),indent=2))
