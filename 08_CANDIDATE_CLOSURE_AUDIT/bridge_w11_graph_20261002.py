from pathlib import Path
import pysam,json,collections,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_GRAPH_BRIDGE_20261002';O.mkdir(exist_ok=False);I=A/'SV_PBGV000012/ROS_expanded_junction_20261001/reads.fastq'
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with pysam.FastaFile(str(A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa')) as f:j=f.fetch('junction_L4707939_R33613308').upper()
with pysam.FastaFile(str(A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa')) as f:flank=rc(f.fetch('NC_051836.1',35342033,35342133).upper())
reads=[]
with pysam.FastxFile(str(I)) as f:
 for r in f:reads.append((r.name.rsplit('_mate',1)[0],r.sequence.upper(),[ord(x)-33 for x in r.quality]))
results=[]
for k in [31,41,51]:
 edges=collections.defaultdict(set)
 for name,seq,q in reads:
  for s,qq in [(seq,q),(rc(seq),q[::-1])]:
   for i in range(len(s)-k):
    mer=s[i:i+k+1]
    if set(mer)<=set('ACGT') and min(qq[i:i+k+1])>=20:edges[mer].add(name)
 left=k//2;start=j[1000-left:1000+(k-left)];target='A'*left+flank[:k-left]
 queue=collections.deque([start]);distance={start:0};parent={};found=False
 while queue:
  node=queue.popleft()
  if node==target:found=True;break
  if distance[node]>=1000:continue
  for base in 'ACGT':
   if len(edges.get(node+base,()))<2:continue
   nxt=node[1:]+base
   if nxt not in distance:distance[nxt]=distance[node]+1;parent[nxt]=(node,base);queue.append(nxt)
  if len(distance)>1000000:raise RuntimeError('state cap')
 path=None
 if found:
  n=target;ext=[]
  while n!=start:n,b=parent[n];ext.append(b)
  path=start+''.join(ext[::-1])
 results.append({'k':k,'reachable':found,'visited_nodes':len(distance),'max_steps':1000,'shortest_graph_steps':distance.get(target),'path':path,'start':start,'target':target})
 del edges
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':results,'limits':['Graph reachability only, not haplotype assembly; repeat switching possible','Homopolymer cycles cannot determine tail length','BQ20 every edge base and>=2read names per edge','Local retained reads only; negative result is limited to graph criteria','No assembler found in checked PATH/conda locations']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'script.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest());print(json.dumps(results))