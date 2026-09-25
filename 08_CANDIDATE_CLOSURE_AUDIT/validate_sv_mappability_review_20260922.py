from pathlib import Path
import csv,json,statistics,hashlib,datetime
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked/coverage_mappability_20260922'
rows=list(csv.DictReader((A/'windows.tsv').open(),delimiter='\t')); lookup={int(r['start0']):r for r in rows};matches=json.loads((A/'matched_windows.json').read_text());summary=json.loads((A/'summary.json').read_text())
assert len(rows)==len(lookup)==400
assert sum(int(r['canonical_probes']) for r in rows)==39997
assert len(matches)==len({r['start0'] for r in matches})==281
ratios=[]
for m in matches:
 r=lookup[m['start0']]; controls=[lookup[s] for s in m['controls_start0']]
 assert r['relation_to_catalogued_span']=='inside' and 3<=len(controls)<=5
 assert len(set(m['controls_start0']))==len(controls)
 for s in controls:
  assert s['relation_to_catalogued_span']=='outside'
  assert abs(float(s['GC_fraction_canonical'])-float(r['GC_fraction_canonical']))<=.01
  assert abs(float(s['sampled_origin_MAPQ30_fraction'])-float(r['sampled_origin_MAPQ30_fraction']))<=.05
 ratio=float(r['mean_read_base_depth'])/statistics.median(float(s['mean_read_base_depth']) for s in controls)
 assert abs(ratio-m['ratio'])<1e-12;ratios.append(ratio)
assert abs(statistics.median(ratios)-summary['matched_ratio_median'])<1e-12
unmatched=sorted(set(s for s,r in lookup.items() if r['relation_to_catalogued_span']=='inside')-{m['start0'] for m in matches});assert len(unmatched)==7
receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='passed_output_consistency_and_matching_checks',unmatched_inside_starts0=unmatched,sha256={n:hashlib.sha256((A/n).read_bytes()).hexdigest() for n in ['windows.tsv','matched_windows.json','summary.json']})
(A/'validation.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
