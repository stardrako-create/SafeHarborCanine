from pathlib import Path
import json,pysam,subprocess,collections,hashlib,datetime,traceback
A=Path(__file__).resolve().parent/'SV_PBGV000010/alignment_chunked';old=A/'junction_competition_20260922_v2';O=A/'joint_alternative_panel_20260925';O.mkdir(exist_ok=False)
R='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa';B='/home/stardrako/miniforge3/envs/atac/bin/bwa'
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
status={'state':'running','started_utc':now()}
def save():(O/'status.json').write_text(json.dumps(status,indent=2))
save()
try:
 intervals=collections.defaultdict(list)
 with pysam.FastaFile(R) as ref:
  for file in [old/'full_reference.sam',A/'junction_review_20260922/clips_full_reference.sam']:
   with pysam.AlignmentFile(str(file)) as sam:
    for r in sam:
     if not r.is_unmapped:intervals[r.reference_name].append((max(0,r.reference_start-1000),min(ref.get_reference_length(r.reference_name),r.reference_end+1000)))
  intervals['NC_049253.1'] += [(4706938,4708938),(33612308,33614308)]
  merged=[]
  for chrom,spans in sorted(intervals.items()):
   out=[]
   for s,e in sorted(spans):
    if out and s<=out[-1][1]:out[-1][1]=max(out[-1][1],e)
    else:out.append([s,e])
   merged.extend((chrom,s,e) for s,e in out)
  meta={}
  with (O/'panel.fa').open('w') as f:
   for i,(chrom,s,e) in enumerate(merged):
    key=f'genomic_{i}';seq=ref.fetch(chrom,s,e).upper();meta[key]={'chrom':chrom,'start0':s,'end0':e};f.write('>'+key+'\n'+seq+'\n')
   with pysam.FastxFile(str(old/'targets.fa')) as targets:
    for r in targets:
     if r.name.startswith('junction_'):f.write('>'+r.name+'\n'+r.sequence.upper()+'\n')
 (O/'genomic_coordinates.json').write_text(json.dumps(meta,indent=2));(O/'reads.fa').write_bytes((old/'reads.fa').read_bytes())
 commands=[[B,'index',str(O/'panel.fa')],[B,'mem','-a','-t','2',str(O/'panel.fa'),str(O/'reads.fa')]]
 (O/'commands.json').write_text(json.dumps(commands,indent=2))
 for cmd,name in zip(commands,['index.stdout','joint.sam']):
  with (O/name).open('w') as out,(O/(name+'.log')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
 by=collections.defaultdict(list)
 with pysam.AlignmentFile(str(O/'joint.sam')) as sam:
  for r in sam:
   if r.is_unmapped:continue
   junction=r.reference_name.startswith('junction_');left=sum(max(0,min(e,1000)-s) for s,e in r.get_blocks() if s<1000);right=sum(max(0,e-max(s,1000)) for s,e in r.get_blocks() if e>1000)
   by[r.query_name].append({'target':r.reference_name,'start0':r.reference_start,'end0':r.reference_end,'AS':r.get_tag('AS'),'NM':r.get_tag('NM'),'cigar':r.cigarstring,'mapq_panel_only':r.mapping_quality,'reverse':r.is_reverse,'junction_spanning':junction and min(left,right)>=20,'genomic':not junction})
 report=[]
 for name,hits in sorted(by.items()):
  j=[h for h in hits if h['junction_spanning']];g=[h for h in hits if h['genomic']];bj=max((h['AS'] for h in j),default=None);bg=max((h['AS'] for h in g),default=None)
  report.append({'read':name,'molecule':name.rsplit('_mate',1)[0],'best_spanning_AS':bj,'best_genomic_AS':bg,'score_margin':None if bj is None or bg is None else bj-bg,'best_spanning_hits':[h for h in j if h['AS']==bj],'best_genomic_hits':[h for h in g if h['AS']==bg]})
 groups=collections.defaultdict(list)
 for r in report:
  if r['best_spanning_AS'] is not None:groups[r['molecule']].append(r['read'])
 result={'utc':now(),'genomic_segments':len(meta),'genomic_bp':sum(v['end0']-v['start0'] for v in meta.values()),'reads_reported':len(report),'spanning_reads':sum(r['best_spanning_AS'] is not None for r in report),'spanning_molecule_names':dict(groups),'read_results':report,'limits':['Joint index contains reported genomic alternatives from prior full-read and clipped-read alignments plus flanks and nine synthetic junctions; not the entire genome','Independent single-read alignments; molecule names group mates but do not establish paired likelihood or unique molecules','Panel MAPQ is not whole-genome MAPQ; score advantage does not establish genotype','BWA remains heuristic; selection of boundary-clipped reads is not an unbiased SV test']}
 (O/'review.json').write_text(json.dumps(result,indent=2));(O/'all_hits.json').write_text(json.dumps(by,indent=2));status.update(state='completed',genomic_segments=len(meta),spanning_reads=result['spanning_reads'],spanning_molecule_names=len(groups),review_sha256=hashlib.sha256((O/'review.json').read_bytes()).hexdigest())
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
status['updated_utc']=now();save();print(json.dumps(status))
