from pathlib import Path
import pysam,json,collections,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_PATH_AUDIT_20261002';O.mkdir(exist_ok=False);I=A/'SV_PBGV000012/ROS_expanded_junction_20261001/reads.fastq';P=A/'W11_GRAPH_BRIDGE_20261002/review.json';path=json.loads(P.read_text())['results'][0]['path'];rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results={}
for n in [32,42,52,75,100,150]:
 targets={path[i:i+n] for i in range(len(path)-n+1)};support=collections.defaultdict(set)
 with pysam.FastxFile(str(I)) as f:
  for r in f:
   seq=r.sequence.upper();q=[ord(c)-33 for c in r.quality];name=r.name.rsplit('_mate',1)[0]
   for s,qq in [(seq,q),(rc(seq),q[::-1])]:
    for i in range(len(s)-n+1):
     mer=s[i:i+n]
     if mer in targets and min(qq[i:i+n])>=20:support[mer].add(name)
 counts=[len(support[path[i:i+n]]) for i in range(len(path)-n+1)]
 results[str(n)]={'window_count':len(counts),'zero_support_positions0':[i for i,c in enumerate(counts) if c==0],'less_than_two_positions0':[i for i,c in enumerate(counts) if c<2],'counts':counts}
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'path_length':len(path),'results':results,'limits':['Exact windows BQ20 throughout, distinct names not molecules','Zero support does not exclude mismatched or low-quality reads','Window support alone does not phase full path or resolve repeat origin','Path length is not insertion size']};(O/'review.json').write_text(json.dumps(out,indent=2));(O/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),I,P]},indent=2));print(json.dumps({n:{'windows':r['window_count'],'zero':len(r['zero_support_positions0']),'less2':len(r['less_than_two_positions0'])} for n,r in results.items()}))