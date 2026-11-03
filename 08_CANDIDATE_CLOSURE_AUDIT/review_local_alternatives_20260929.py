from pathlib import Path
import json,pysam,datetime,hashlib,collections
A=Path(__file__).resolve().parent;S=A/'SV_CALL_AUDIT_20260926/nearby_nonreference_calls.json';O=A/'W11_LOCAL_ALTERNATIVES_20260929';O.mkdir(exist_ok=False)
raw=json.loads(S.read_text());groups=collections.defaultdict(list)
for r in raw['rows']:groups[r['id'].split('.')[0]].append(r)
B=A/'SV_PBGV000010/alignment_chunked/selected_markdup.bam';reports=[]
with pysam.AlignmentFile(str(B),'rb') as bam:
 for key,models in groups.items():
  r=next((x for x in models if x['info'].get('SVMODEL')=='AGGREGATED'),models[0]);pos=r['pos'];events=[];names=set();span=set()
  for read in bam.fetch('NC_049253.1',max(0,pos-201),pos+200):
   if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail or read.mapping_quality<30:continue
   names.add(read.query_name)
   if any(s<=pos-20 and e>=pos+20 for s,e in read.get_blocks()):span.add(read.query_name)
   cursor=read.reference_start
   for i,(op,n) in enumerate(read.cigartuples or []):
    if op in (1,2,4) and n>=20 and abs(cursor-pos)<=200:
     events.append({'read':read.query_name,'op':{1:'I',2:'D',4:'S'}[op],'pos0':cursor,'length':n,'mapq':read.mapping_quality,'cigar':read.cigarstring})
    if op in (0,2,3,7,8):cursor+=n
  reports.append({'event_id':key,'pos1':pos,'end1':r['end'],'alts':r['alts'],'models':[{'id':x['id'],'model':x['info'].get('SVMODEL'),'site_filters':x['filters'],'sample_FT':x['sample'].get('FT'),'GT':x['sample'].get('GT')} for x in models],'distance_to_nearest_catalogue_boundary_bp':min(abs(pos-4707939),abs(pos-33613308)),'eligible_read_names':len(names),'reference_block_spanning_names':len(span),'cigar_events':events})
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(S.read_bytes()).hexdigest(),'record_count':len(raw['rows']),'event_groups':reports,'method':'Group model suffixes by base record ID. Inspect primary nonduplicate QC-passing MAPQ>=30 alignments within200bp of each POS. Report I/D/S>=20bp near POS; descriptive screen, not variant calling.','limits':['Grouping is catalogue identity, not demonstrated biological independence','Missing FT is unknown, not PASS','CIGAR events are not necessarily support for symbolic ALT; small/complex events and low MAPQ reads excluded','Regional retention and read names do not establish independent molecules','Calls several kb from a boundary do not directly explain the seven selected 150bp reads; no causal inference']}
(O/'review.json').write_text(json.dumps(res,indent=2))
for x in reports:print(x['event_id'],'distance',x['distance_to_nearest_catalogue_boundary_bp'],'models',[(m['model'],m['sample_FT']) for m in x['models']],'names',x['eligible_read_names'],'events',[(e['op'],e['pos0'],e['length']) for e in x['cigar_events']])
