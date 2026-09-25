from pathlib import Path
import json,hashlib,collections,datetime,pysam
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=C/'pending_allele_diagnostics_20260923';O.mkdir(exist_ok=False)
integrated=json.loads((C/'integrated_review_20260923/review.json').read_text());checks={f:hashlib.sha256((C/f).read_bytes()).hexdigest()==h for f,h in integrated['sources_sha256'].items()};assert all(checks.values())
base=json.loads((C/'gapped_projection_20260923/UU_base_maps.json').read_text());windows={r['id']:r for r in json.loads((C/'gapped_projection_20260923/UU_review.json').read_text())['results']};pending=json.loads((C/'partial_variant_projection_20260923_v3/pending.json').read_text());rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results=[]
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'
with pysam.FastaFile(str(R)) as ref:
 for item in pending:
  a=item['allele'];key=item['id'].split('__')[0];w=windows[key];result={'id':item['id'],'source_allele':a}
  if not a['sequence_allele']:result['diagnosis']='nonsequence_ALT';results.append(result);continue
  inv=collections.defaultdict(list)
  for i,v in enumerate(base[key]):
   if v is not None:inv[v['uu_pos0']].append(i)
  coords=list(range(a['pos1']-1,a['pos1']-1+len(a['ref'])))
  if any(len(inv[p])!=1 for p in coords):result['diagnosis']='source_REF_contains_unmapped_or_nonunique_base';results.append(result);continue
  ix=[inv[p][0] for p in coords];step=1 if w['strand']=='+' else -1
  if any(ix[i]!=ix[0]+i*step for i in range(len(ix))):result['diagnosis']='source_REF_crosses_ROS_insertion_or_discontinuity';results.append(result);continue
  lo=w['ros_start0']+min(ix);hi=w['ros_start0']+max(ix)+1;ros=ref.fetch(w['ros_chrom'],lo,hi).upper();uuref=a['ref'].upper();alt=a['alt'].upper()
  if w['strand']=='-':uuref=rc(uuref);alt=rc(alt)
  assert ros!=uuref,'Unexpected exactly matching REF in pending set'
  result.update(ros_chrom=w['ros_chrom'],ros_start0=lo,ros_end0=hi,ros_reference=ros,oriented_source_ref=uuref,oriented_source_alt=alt)
  result['diagnosis']='source_ALT_equals_ROS_reference' if alt==ros else 'REF_disagreement_ALT_also_differs_from_ROS'
  results.append(result)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_pending_records':len(pending),'diagnosis_counts':dict(collections.Counter(r['diagnosis'] for r in results)),'integrated_source_hash_checks':checks,'results':results,'limitations':['Diagnostic aligned-coordinate comparison only; no previously rejected variant promoted','Source ALT frequency is not converted into ROS ALT frequency or complemented','ALT matching ROS reference may reflect reference-allele differences; does not establish donor genotype or absence of variation','Indel/phase reconciliation remains unresolved']}
(O/'review.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary['diagnosis_counts'],indent=2))
