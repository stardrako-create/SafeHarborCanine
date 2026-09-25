from pathlib import Path
import pysam,json,collections,datetime
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked';O=A/'local_extension_20260922';O.mkdir(exist_ok=True)
ref='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa';chrom='NC_049253.1';K=51
comp=str.maketrans('ACGTN','TGCAN');rc=lambda s:s.translate(comp)[::-1]
reads={}
with pysam.AlignmentFile(str(A/'selected_markdup.bam'),'rb') as f:
 for pos in [4707939,33613308]:
  for r in f.fetch(chrom,pos-2000,pos+2000):
   if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality<20:continue
   key=r.query_name+('_1' if r.is_read1 else '_2');reads[key]=(r.query_name,r.query_sequence,list(r.query_qualities))
edges=collections.defaultdict(set)
for key,(name,seq,q) in reads.items():
 for s,qual in [(seq,q),(rc(seq),q[::-1])]:
  for i in range(len(s)-K):
   frag=s[i:i+K+1]
   if 'N' not in frag and min(qual[i:i+K+1])>=20:edges[frag].add(name)
with pysam.FastaFile(ref) as f:
 seed=f.fetch(chrom,4707939-K,4707939);left=f.fetch(chrom,4707939,4708939);right=f.fetch(chrom,33613308,33614308)
 # Seed has a fixed reference coordinate; counts below are graph support, not phased molecules.
 graph={}
 for k,names in edges.items():
  if len(names)>=2:graph.setdefault(k[:-1],[]).append((k[-1],len(names)))
 paths=[('',seed,[],{seed})];finished=[];states=0;maxstates=20000
 while paths and states<maxstates:
  seq,node,support,visited=paths.pop();states+=1
  options=graph.get(node,[])
  if len(seq)>=500 or not options:
   finished.append(dict(extension=seq,edge_support=support,stop='length_limit' if len(seq)>=500 else 'no_supported_edge'));continue
  for base,n in options:
   nxt=node[1:]+base
   if nxt in visited:finished.append(dict(extension=seq+base,edge_support=support+[n],stop='repeat_cycle'));continue
   paths.append((seq+base,nxt,support+[n],visited|{nxt}))
 def prefix(a,b):
  for i,(x,y) in enumerate(zip(a,b)):
   if x!=y:return i
  return min(len(a),len(b))
 for x in finished:
  x['length']=len(x['extension']);x['match_left_reference_prefix']=prefix(x['extension'],left);x['match_distal_reference_prefix']=prefix(x['extension'],right)
  x['minimum_edge_support']=min(x['edge_support']) if x['edge_support'] else 0
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),method='Exploratory exact 51-mer graph; each edge requires >=2 unique read names and Q>=20 across 52 bases; MAPQ>=20 nonduplicate primary local reads',reads=len(reads),unique_names=len({v[0] for v in reads.values()}),seed=seed,seed_coordinate_end0=4707939,states=states,unfinished_paths=len(paths),paths=finished,limitations=['Graph paths can combine different molecules/haplotypes; not a validated assembly or SV genotype','Reads retained near two preselected loci; not genome-wide local assembly','Repeat cycles terminate paths; absence of extension is not absence of variant','Quality, length and minimum edge-support filters are explicit exploratory choices'])
(O/'graph_extensions.json').write_text(json.dumps(result,indent=2));(O/'extensions.fa').write_text(''.join('>path'+str(i)+'\n'+seed+x['extension']+'\n' for i,x in enumerate(finished)))
print(json.dumps({k:v for k,v in result.items() if k!='paths'},indent=2));print(json.dumps([{k:v for k,v in x.items() if k not in ['extension','edge_support']} for x in finished],indent=2))
