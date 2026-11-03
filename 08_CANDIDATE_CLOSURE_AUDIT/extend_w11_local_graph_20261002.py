from pathlib import Path
import pysam,json,collections,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_LOCAL_EXTENSION_20261002';O.mkdir(exist_ok=False)
I=A/'SV_PBGV000012/ROS_expanded_junction_20261001/reads.fastq'
J=A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa'
with pysam.FastaFile(str(J)) as f:j=f.fetch('junction_L4707939_R33613308').upper()
reads=[]
with pysam.FastxFile(str(I)) as f:
 for r in f:reads.append((r.name.rsplit('_mate',1)[0],r.sequence.upper(),[ord(x)-33 for x in r.quality]))
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results=[]
for k in (31,41,51):
 edges=collections.defaultdict(set)
 for name,seq,qual in reads:
  for ss,qq in ((seq,qual),(rc(seq),qual[::-1])):
   for i in range(len(ss)-k):
    mer=ss[i:i+k+1]
    if set(mer)<=set('ACGT') and min(qq[i:i+k+1])>=20:edges[mer].add(name)
 # A sequence spanning the already observed junction, with 40bp each side.
 seed=j[960:1040];walk=seed;trace=[];seen=set();reason='length_cap'
 for step in range(300):
  suffix=walk[-k:]
  if suffix in seen:reason='cycle';break
  seen.add(suffix)
  counts={b:len(edges.get(suffix+b,())) for b in 'ACGT'}
  supported=[b for b,n in counts.items() if n>=2]
  trace.append({'right_position_next0':len(walk)-40,'counts_distinct_names':counts})
  if len(supported)!=1:reason='branch' if len(supported)>1 else 'insufficient_support';break
  walk+=supported[0]
 results.append({'k':k,'seed_left_bp':40,'seed_right_bp':40,'extended_right_bp':len(walk)-40,'stop':reason,'sequence':walk,'trace':trace,'exact_model_prefix_bp':next((i for i,(x,y) in enumerate(zip(walk[40:],j[1000:])) if x!=y),len(walk)-40)})
 del edges
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_reads':len(reads),'results':results,'limits':['Exploratory local de Bruijn traversal, not resolved assembly','Requires at least2 distinct read names per edge and BQ20 on every base','Read names are not independently confirmed molecules','Edges may switch between repeat copies; extension is not proof of physical linkage','Seed contains40bp each side of observed junction; local subset and retained mates only','Branches stop traversal rather than selecting highest count']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),I]},indent=2));print(json.dumps([{k:v for k,v in x.items() if k not in ['trace','sequence']} for x in results]))