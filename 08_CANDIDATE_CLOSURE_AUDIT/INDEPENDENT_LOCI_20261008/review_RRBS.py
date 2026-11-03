from pathlib import Path
import json,pyBigWig,pysam,math
O=Path(__file__).parent;P=O.parent.parent
rows=json.loads((O/'TAD_clear_queue.json').read_text(encoding='utf-8-sig'));out=[]
with pysam.FastaFile(str(P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa')) as fa,pyBigWig.open(str(P/'06_v1.20.0/RRBS/methylation_weighted_mean.bw')) as bw,pyBigWig.open(str(P/'06_v1.20.0/RRBS/cpg_coverage_frequency.bw')) as cov:
 for w in rows:
  c=w['chrom'];s=int(w['start0']);e=int(w['end0']);seq=fa.fetch(c,s,e)
  v=bw.values(c,s,e);d=cov.values(c,s,e);observed=[x for x,n in zip(v,d) if math.isfinite(n) and n>0 and math.isfinite(x)]
  out.append(dict(w,softmasked_bp=sum(x.islower() for x in seq),RRBS_observed_bp=len(observed),RRBS_observed_mean=sum(observed)/len(observed) if observed else None))
(O/'sequence_RRBS_review.json').write_text(json.dumps(dict(results=out,limitations=['Softmask is reference annotation, not fresh RepeatMasker','Coverage-positive finite observations only; uncovered bases not zero methylation']),indent=2))
print(json.dumps(sorted(out,key=lambda x:(x['softmasked_bp'], -float(x['atac_percentile'])))[:5]))