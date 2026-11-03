from pathlib import Path
import json,itertools,pysam,hashlib
N=Path(__file__).resolve().parent;S=N/'hypothetical_TTTV23_panel_20261009/results.json';I=N/'ANO2_indel_union_robustness_20261009.json';w='NC_051831.1_39708074_39709074';start=39708074
ref=pysam.FastaFile(str(N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa')).fetch('NC_051831.1',start,start+1000).upper();indels=json.loads(I.read_text())['indels'];guides=[g for g in json.loads(S.read_text())['results'] if g['window']==w];assert len(guides)==14
checks=[]
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
 dna.extend(ref[cursor:]);origins.extend(range(cursor,1000));dna=''.join(dna);assert len(dna)==len(origins);lookup={p:i for i,p in enumerate(origins) if p is not None};failed=[]
 for g in guides:
  off=g['start0']-start;lo=off-4 if g['strand']=='+' else off;j=lookup.get(lo)
  okay=j is not None and origins[j:j+27]==list(range(lo,lo+27)) and dna[j:j+27]==ref[lo:lo+27]
  if not okay:failed.append(dict(start0=g['start0'],strand=g['strand']))
 checks.append(dict(selected_indels=[r['id'] for use,r in zip(selected,indels) if use],guides_checked=len(guides),failed=failed,alternate_sequence_sha256=hashlib.sha256(dna.encode()).hexdigest()))
assert all(not r['failed'] for r in checks)
out=N/'ANO2_hypothetical_TTTV23_PAM_indel_robustness_20261009.json';assert not out.exists();out.write_text(json.dumps(dict(state='completed',window=w,guides_with_intact_protospacer_and_PAM_all32_combinations=14,checks=checks,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S,I]},limits=['Synthetic union of five PASS local indels, not observed haplotypes or donor genotypes.','Canonical TTTV hypothetical panel; no guide approval or alternate-genome off-target clearance.','Local SNPs/REF discrepancies masked before panel construction.']),indent=2));print(json.dumps(dict(guides=14,combinations=32,failed=0)))

