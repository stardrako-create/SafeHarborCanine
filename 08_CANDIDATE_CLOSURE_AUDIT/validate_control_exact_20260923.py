from pathlib import Path
import csv,json,hashlib,collections,datetime,pysam
O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922';D=O/'exact_sequences';r=json.loads((D/'exact_sequence_review.json').read_text());state=json.loads((D/'exact_scan_status.json').read_text());rows=list(csv.DictReader((D/'all_short_segments_exact_uniqueness.tsv').open(),delimiter='\t'))
assert r['status']==state['status']=='completed'
with pysam.FastaFile(r['reference']) as f:
 assert sum(f.lengths)==r['bases_scanned']==state['bases_scanned'];assert len(f.references)==r['contigs']
assert len(rows)==76180;assert len({(x['window_id'],x['start0'],x['length']) for x in rows})==len(rows)
groups=collections.defaultdict(list)
for x in rows:
 assert int(x['exact_locus_count_both_strands'])>=1;assert (x['unique_exact']=='True')==(int(x['exact_locus_count_both_strands'])==1);groups[x['window_id'],int(x['length'])].append(x)
assert len(groups)==80
for x in r['summary']:
 v=groups[x['window_id'],x['length']];assert len(v)==1001-x['length']==x['segments'];assert sum(y['unique_exact']=='True' for y in v)==x['unique'];assert len(v)-x['unique']==x['nonunique'];assert max(int(y['exact_locus_count_both_strands']) for y in v)==x['max_locus_count']
tad=json.loads((O/'TAD_multiscale_review.json').read_text());epic=json.loads((O/'EpiC_roundtrip_and_state_review.json').read_text());eligible=set(tad['no_hit_all_scales_ids'])&{x['id'] for x in epic['state_review'] if x['roundtrip_pass'] and not x['promoter_enhancer_tissues']}
summary=[x for x in r['summary'] if x['window_id'] in eligible and x['length']==20]
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed_internal_consistency_and_reference_extent_checks',segments=len(rows),bases_scanned=r['bases_scanned'],eligible_after_TAD_EpiC_ids=sorted(eligible),twenty_base_summary=summary,sha256={n:hashlib.sha256((D/n).read_bytes()).hexdigest() for n in ['exact_sequence_review.json','all_short_segments_exact_uniqueness.tsv']},limitations=['Consistency verification and existing brute-force toy self-test; not an independent full-genome reimplementation','Exact counts do not establish PAM suitability, variant clearance or mismatch/bulge specificity'])
(D/'validation_20260923.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
