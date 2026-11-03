from pathlib import Path
import json,pysam,itertools,hashlib
N=Path(__file__).resolve().parent;S=N/'normalized_complete_20261008/review.json';rows=json.loads(S.read_text())['results'];F=N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';fa=pysam.FastaFile(str(F));panels={name:json.loads((N/path).read_text())['results'] for name,path in [('NGG20','hypothetical_NGG_panel_20261009/results.json'),('TTTV23','hypothetical_TTTV23_panel_20261009/results.json')]};results=[]
for w in ['NC_051807.1_2091307_2092307','NC_051814.1_43352999_43353999']:
 chrom,start,end=w.rsplit('_',2);start=int(start);ref=fa.fetch(chrom,start,int(end)).upper();indels=sorted([r for r in rows if r['window']==w and r['source_allele']['filters']==['PASS'] and len(r['ref'])!=len(r['alt'])],key=lambda r:r['start0']);assert len(indels)<10
 for r in indels:assert r['haplotype_equivalence_verified'] and ref[r['start0']-start:r['start0']-start+len(r['ref'])]==r['ref'].upper()
 checks=[]
 for bits in itertools.product([False,True],repeat=len(indels)):
  dna=[];origins=[];cursor=0
  for use,r in zip(bits,indels):
   lo=r['start0']-start;hi=lo+len(r['ref']);dna.extend(ref[cursor:lo]);origins.extend(range(cursor,lo))
   if use:
    a=r['ref'].upper();b=r['alt'].upper();prefix=0
    while prefix<min(len(a),len(b)) and a[prefix]==b[prefix]:prefix+=1
    suffix=0
    while suffix<min(len(a),len(b))-prefix and a[-1-suffix]==b[-1-suffix]:suffix+=1
    dna.extend(b);origins.extend(list(range(lo,lo+prefix))+[None]*(len(b)-prefix-suffix)+list(range(hi-suffix,hi)))
   else:dna.extend(ref[lo:hi]);origins.extend(range(lo,hi))
   cursor=hi
  dna.extend(ref[cursor:]);origins.extend(range(cursor,1000));dna=''.join(dna);lookup={p:i for i,p in enumerate(origins) if p is not None};tested={};failed=[]
  for name,rr in panels.items():
   gg=[g for g in rr if g['window']==w];tested[name]=len(gg)
   for g in gg:
    off=g['start0']-start
    if name=='NGG20':lo=off if g['strand']=='+' else off-3;length=23
    else:lo=off-4 if g['strand']=='+' else off;length=27
    j=lookup.get(lo);okay=j is not None and origins[j:j+length]==list(range(lo,lo+length)) and dna[j:j+length]==ref[lo:lo+length]
    if not okay:failed.append(dict(panel=name,start0=g['start0'],strand=g['strand']))
  checks.append(dict(indels_selected=[r['id'] for use,r in zip(bits,indels) if use],guides_tested=tested,failed=failed,sequence_sha256=hashlib.sha256(dna.encode()).hexdigest()))
 assert all(not c['failed'] for c in checks);results.append(dict(window=w,PASS_indels=indels,synthetic_combinations=len(checks),checks=checks))
out=N/'independent_two_loci_hypothetical_PAM_indel_robustness_20261009.json';assert not out.exists();out.write_text(json.dumps(dict(state='completed',results=results,source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),limits=['Observed local PASS indels only, SNP/REF masks applied before constructing panels; no donor validation or populationwidealternate-genomeofftargets.','No indel observation is not evidence of no individual indels;1900 has zero known PASS indels in this query.','HypotheticalPAMpanels not actual nuclease or guide approval.']),indent=2));print(json.dumps([dict(window=r['window'],indels=len(r['PASS_indels']),combinations=r['synthetic_combinations'],tested=r['checks'][0]['guides_tested']) for r in results]))
