from pathlib import Path
import json,csv,pysam,itertools,hashlib,collections
N=Path(__file__).resolve().parent;A=N.parent/'ALL_CANDIDATES_20261005';P=N.parent.parent;w='NC_051831.1_39708074_39709074';chrom,start,end=w.rsplit('_',2);start=int(start);end=int(end)
R=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';fa=pysam.FastaFile(str(R));ref=fa.fetch(chrom,start,end).upper();mask=set(next(r for r in json.loads((A/'masked_specificity_20261008.json').read_text())['results'] if r['id']==w)['masked_offsets'])
source=A/'normalized_exact_blocks_20261006_v2/review.json';indels=[r for r in json.loads(source.read_text())['results'] if r['window']==w and r['source_allele']['filters']==['PASS'] and len(r['ref'])!=len(r['alt'])];indels.sort(key=lambda r:r['start0']);assert len(indels)==5
for a,b in zip(indels,indels[1:]):assert a['start0']+len(a['ref'])<=b['start0']
for r in indels:assert r['haplotype_equivalence_verified'] and ref[r['start0']-start:r['start0']-start+len(r['ref'])]==r['ref'].upper()
segments=[];tsv=A/'exact_specificity_20261006/all_short_segments_exact_uniqueness.tsv'
for r in csv.DictReader(tsv.open(),delimiter='\t'):
 if r['window_id']==w and int(r['length'])==20 and r['unique_exact']=='True':
  off=int(r['start0'])-start
  if not any(i in mask for i in range(off,off+20)):segments.append(off)
assert len(segments)==631
results=[]
for selected in itertools.product([False,True],repeat=5):
 dna=[];origins=[];cursor=0
 for use,r in zip(selected,indels):
  lo=r['start0']-start;hi=lo+len(r['ref']);dna.extend(ref[cursor:lo]);origins.extend(range(cursor,lo))
  if use:
   a=r['ref'].upper();b=r['alt'].upper();prefix=0
   while prefix<min(len(a),len(b)) and a[prefix]==b[prefix]:prefix+=1
   suffix=0
   while suffix<min(len(a),len(b))-prefix and a[-1-suffix]==b[-1-suffix]:suffix+=1
   dna.extend(b);origins.extend(list(range(lo,lo+prefix))+[None]*(len(b)-prefix-suffix)+list(range(hi-suffix,hi)))
  else:dna.extend(ref[lo:hi]);origins.extend(range(lo,hi))
  cursor=hi
 dna.extend(ref[cursor:]);origins.extend(range(cursor,1000));dna=''.join(dna);assert len(dna)==len(origins);lookup={p:i for i,p in enumerate(origins) if p is not None};passed=[];failed=[]
 for off in segments:
  idx=lookup.get(off)
  okay=idx is not None and origins[idx:idx+20]==list(range(off,off+20)) and dna[idx:idx+20]==ref[off:off+20]
  (passed if okay else failed).append(off)
 results.append(dict(selected_indel_ids=[r['id'] for use,r in zip(selected,indels) if use],length=len(dna),sequence_sha256=hashlib.sha256(dna.encode()).hexdigest(),eligible_segments_intact=len(passed),failed_offsets=failed))
assert all(not r['failed_offsets'] for r in results)
out=N/'ANO2_indel_union_robustness_20261009.json';assert not out.exists();out.write_text(json.dumps(dict(window=w,indels=indels,eligible_20bp_segments=631,synthetic_combinations=32,results=results,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),limitations=['Synthetic combinations of five non-overlapping PASS indels, not observed/phased haplotypes.','Checks unchanged sequence at its intended site outside existing PASS REF-span mask; not genome-wide alternate-allele off-targets or PAM/bulges.','Symbolic star associations retained as unphased; union includes their candidate deletions without assuming cis.','No biological or editing-safety approval.']),indent=2));print(json.dumps(dict(indels=5,synthetic_combinations=32,segments_preserved_in_every_combination=631)))