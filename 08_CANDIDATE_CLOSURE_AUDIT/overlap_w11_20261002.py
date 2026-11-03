from pathlib import Path
import pysam,json,datetime,hashlib,collections
A=Path(__file__).resolve().parent;I=A/'SV_PBGV000012/ROS_expanded_junction_20261001';O=A/'W11_OVERLAP_20261002';O.mkdir(exist_ok=False)
review=json.loads((I/'review.json').read_text());names={r['read'].rsplit('_mate',1)[0] for r in review['rows']};keys={r['read'] for r in review['rows']}
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
reads=[]
with pysam.FastxFile(str(I/'reads.fastq')) as f:
 for r in f:
  q=[ord(c)-33 for c in r.quality]
  reads.extend([(r.name,'+',r.sequence.upper(),q),(r.name,'-',rc(r.sequence.upper()),q[::-1])])
seeds=[]
with pysam.AlignmentFile(str(I/'joint.sam'),'r') as f:
 for r in f:
  if r.query_name in keys and r.reference_name=='junction_L4707939_R33613308':
   seeds.append((r.query_name,r.query_sequence,r.query_qualities,r.reference_end))
rows=[]
for name,s,q,end in seeds:
 candidates=[]
 for other,strand,t,qual in reads:
  if other==name:continue
  for n in range(min(len(s),len(t))-1,59,-1):
   if s[-n:]==t[:n] and min(q[-n:])>=20 and min(qual)>=20:
    candidates.append({'read':other,'strand':strand,'overlap':n,'extension':t[n:],'right_end_if_correct':end-1000+len(t)-n,'same_pair_as_seed':other.rsplit('_mate',1)[0]==name.rsplit('_mate',1)[0],'belongs_to_any_spanning_pair':other.rsplit('_mate',1)[0] in names});break
 rows.append({'seed':name,'seed_right_end':end-1000,'candidates':candidates})
summary={}
for n in [60,80,100]:
 hits=[c for r in rows for c in r['candidates'] if c['overlap']>=n and c['right_end_if_correct']>122]
 summary[str(n)]={'candidate_overlaps_beyond122':len(hits),'distinct_reads':len({x['read'] for x in hits}),'same_pair_links':sum(x['same_pair_as_seed'] for x in hits),'spanning_pair_reads':len({x['read'] for x in hits if x['belongs_to_any_spanning_pair']}),'max_hypothetical_right_end':max([x['right_end_if_correct'] for x in hits],default=None)}
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'summary':summary,'rows':rows,'limits':['Exploratory single-step exact overlaps; not an assembly','BQ20 throughout incoming read and seed overlap; mismatched or lower-quality continuations excluded','Repeated overlaps do not establish physical linkage','Pair membership alone does not establish insert orientation or structural support','Candidates beyond122 are hypotheses, not observed spanning reads']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'provenance.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),I/'reads.fastq',I/'joint.sam']},indent=2));print(json.dumps(summary))