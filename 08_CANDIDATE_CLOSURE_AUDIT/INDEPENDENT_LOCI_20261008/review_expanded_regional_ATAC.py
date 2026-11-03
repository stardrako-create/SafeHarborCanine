from pathlib import Path
import json,sys,hashlib,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent
sys.path[:0]=[str(P/'06_v1.21.8/code'),str(P/'06_v1.21.1/code')]
import score_ship_candidates_v2 as legacy
import pyBigWig
cmd=json.loads((P/'06_v1.21.8/scoring_command.json').read_text());track=Path(cmd[cmd.index('--atac-mean-bw')+1]);threshold=float(cmd[cmd.index('--min-atac-accessibility-percentile')+1]);bw=pyBigWig.open(str(track))
rows=json.loads((N/'flanking_coding_context.json').read_text())['results'];backgrounds={};results=[]
for r in rows:
 mid=(r['window_start0']+r['window_end0'])//2
 intervals=[('full_coding_gene_gap',r['intergenic_start0'],r['intergenic_end0']),('centered_50kb_sensitivity',mid-25000,mid+25000),('centered_75kb_sensitivity',mid-37500,mid+37500)]
 for kind,s,e in intervals:
  length=e-s
  if length not in backgrounds:backgrounds[length]=legacy.build_matched_window_background(bw,[length])
  mean=legacy.track_value(bw,r['chrom'],s,e);pct=legacy.percentile_rank(backgrounds[length],mean)
  results.append(dict(id=r['id'],kind=kind,chrom=r['chrom'],start0=s,end0=e,length=length,mean=mean,percentile=pct,above_legacy_accessibility_floor=pct is not None and pct>=threshold,convergent=r['convergent'],original_ship_generation_eligible=kind=='full_coding_gene_gap' and r['convergent'] and 50000<=length<=75000))
out=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),results=results,threshold=threshold,backgrounds={str(k):v.tolist() for k,v in backgrounds.items()},track_sha256=hashlib.sha256(track.read_bytes()).hexdigest(),legacy_code_sha256=hashlib.sha256(Path(legacy.__file__).read_bytes()).hexdigest(),interpretation='Diagnostic for independent expansion, not a substitute for original SHIP regional PASS. Full coding-gene gaps use the same interval definition but differ in size/orientation from original generation; centered windows are sensitivity analyses only. No candidate approval follows from these measurements.')
path=N/'expanded_regional_ATAC_diagnostic_20261009.json';assert not path.exists();path.write_text(json.dumps(out,indent=2));print(json.dumps(results))