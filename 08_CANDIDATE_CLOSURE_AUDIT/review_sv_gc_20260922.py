from pathlib import Path
import csv,json,statistics,datetime
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked'; O=A/'junction_review_20260922'
clips={r['id']:r for r in json.loads((O/'clips.json').read_text())};hits=json.loads((O/'clip_hits.json').read_text()); reps=[]
for cid in sorted({h['clip_id'] for h in hits if h['near_opposite_boundary']}):
 exact=[h for h in hits if h['clip_id']==cid and h['NM']==0 and h['aligned_query_bases']==clips[cid]['clip_length']]
 reps.append(dict(clip_id=cid,read_name=clips[cid]['name'],clip_length=clips[cid]['clip_length'],reported_full_length_exact_alignments=len(exact),distinct_contigs=len({h['chrom'] for h in exact}),anchor_end0=clips[cid]['anchor_end0']))
rows=list(csv.DictReader((A/'coverage_windows_descriptive.tsv').open(),delimiter='\t'))
rows=[r for r in rows if r['MAPQ_min']=='30' and r['exclude_marked_duplicates']=='True' and int(r['bases'])==100000 and float(r['canonical_fraction'])>=.99]
inside=[r for r in rows if r['relation_to_catalogued_span']=='inside'];outside=[r for r in rows if r['relation_to_catalogued_span']=='outside']; matched=[]
for r in inside:
 controls=sorted(outside,key=lambda s:abs(float(s['GC_fraction_canonical'])-float(r['GC_fraction_canonical'])))[:5]
 controls=[s for s in controls if abs(float(s['GC_fraction_canonical'])-float(r['GC_fraction_canonical']))<=.01]
 if len(controls)<3:continue
 base=statistics.median(float(s['mean_read_base_depth']) for s in controls)
 if base>0:matched.append(dict(start0=int(r['start0']),control_count=len(controls),ratio=float(r['mean_read_base_depth'])/base,controls_start0=[int(s['start0']) for s in controls]))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),exact_clip_ambiguity=reps,gc_matched=dict(method='3-5 nearest outside windows within 0.01 absolute GC fraction; controls reusable; MAPQ30 exclude duplicates',inside_eligible=len(inside),inside_matched=len(matched),median_ratio=statistics.median(r['ratio'] for r in matched),windows=matched),limitations=['Reported BWA alignments are not exhaustive genomic exact-match counts','GC matching does not match mappability and outside windows are not proven diploid','Reused controls and correlated windows preclude treating these as independent replicates','No definitive SV genotype assigned'])
(O/'ambiguity_and_gc_review.json').write_text(json.dumps(result,indent=2)); result['gc_matched'].pop('windows');print(json.dumps(result,indent=2))
