"""Sampled reference-remapping proxy; not exhaustive or donor-specific mappability."""
from pathlib import Path
import csv,json,subprocess,pysam,collections,statistics,datetime
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked'; O=A/'coverage_mappability_20260922';O.mkdir(exist_ok=True)
REF='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa';BWA='/home/stardrako/miniforge3/envs/atac/bin/bwa'
def status(stage):
 (O/'status.json').write_text(json.dumps(dict(stage=stage,utc=datetime.datetime.now(datetime.timezone.utc).isoformat())))
status('prepare_reference_probes')
rows=list(csv.DictReader((A/'coverage_windows_descriptive.tsv').open(),delimiter='\t'))
rows=[r for r in rows if r['MAPQ_min']=='30' and r['exclude_marked_duplicates']=='True' and int(r['bases'])==100000 and float(r['canonical_fraction'])>=.99 and r['relation_to_catalogued_span'] in ('inside','outside')]
meta={};counts=collections.defaultdict(collections.Counter)
with pysam.FastaFile(REF) as ref,(O/'probes.fa').open('w') as out:
 for i,r in enumerate(rows):
  for offset in range(0,100000-149,1000):
   start=int(r['start0'])+offset;seq=ref.fetch(r['chrom'],start,start+150).upper();counts[i]['sampled']+=1
   if len(seq)!=150 or set(seq)-set('ACGT'):counts[i]['excluded_noncanonical']+=1;continue
   name=f'w{i}_p{offset}';meta[name]=(i,r['chrom'],start);out.write('>'+name+'\n'+seq+'\n');counts[i]['canonical_probes']+=1
cmd=[BWA,'mem','-t','4',REF,str(O/'probes.fa')];(O/'command.json').write_text(json.dumps(cmd));status('align_reference_probes')
with (O/'probes.sam').open('w') as out,(O/'bwa.stderr.log').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True)
seen=set()
with pysam.AlignmentFile(str(O/'probes.sam'),'r') as f:
 for r in f:
  if r.is_secondary or r.is_supplementary:continue
  assert r.query_name not in seen;seen.add(r.query_name);i,c,s=meta[r.query_name]
  if not r.is_unmapped and r.reference_name==c and r.reference_start==s and r.cigarstring=='150M' and r.mapping_quality>=30:counts[i]['origin_MAPQ30']+=1
assert len(seen)==len(meta)
for i,r in enumerate(rows):
 r.update(counts[i]);r['sampled_origin_MAPQ30_fraction']=counts[i]['origin_MAPQ30']/counts[i]['canonical_probes'] if counts[i]['canonical_probes'] else None
with (O/'windows.tsv').open('w',newline='') as f:
 fields=list(dict.fromkeys(k for r in rows for k in r));w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rows)
inside=[r for r in rows if r['relation_to_catalogued_span']=='inside'];outside=[r for r in rows if r['relation_to_catalogued_span']=='outside'];matches=[]
for r in inside:
 candidates=[s for s in outside if abs(float(s['GC_fraction_canonical'])-float(r['GC_fraction_canonical']))<=.01 and abs(s['sampled_origin_MAPQ30_fraction']-r['sampled_origin_MAPQ30_fraction'])<=.05]
 candidates.sort(key=lambda s:abs(float(s['GC_fraction_canonical'])-float(r['GC_fraction_canonical']))/.01+abs(s['sampled_origin_MAPQ30_fraction']-r['sampled_origin_MAPQ30_fraction'])/.05)
 controls=candidates[:5]
 if len(controls)<3:continue
 depth=statistics.median(float(s['mean_read_base_depth']) for s in controls)
 if depth<=0:continue
 matches.append(dict(start0=int(r['start0']),ratio=float(r['mean_read_base_depth'])/depth,controls_start0=[int(s['start0']) for s in controls]))
(O/'matched_windows.json').write_text(json.dumps(matches,indent=2))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),method='Error-free reference 150-mers at 1kb spacing; fraction remapping exactly to origin at MAPQ30; not exhaustive mappability',eligible_windows=len(rows),canonical_probes=len(meta),inside_eligible=len(inside),inside_matched=len(matches),matched_ratio_median=statistics.median(r['ratio'] for r in matches) if matches else None,control_criteria='3-5 external windows within absolute GC 0.01 and sampled mapping fraction 0.05; controls reused',limitations=['Reference-derived single-end probes omit individual variation, sequencing error and paired-end rescue','Regular 1kb sampling is a proxy, not exhaustive uniqueness','Outside regions not confirmed diploid; controls reused and windows correlated','No formal CNV likelihood or SV genotype; coverage cannot resolve repeat-associated junction'])
(O/'summary.json').write_text(json.dumps(result,indent=2));status('completed');print(json.dumps(result,indent=2))
