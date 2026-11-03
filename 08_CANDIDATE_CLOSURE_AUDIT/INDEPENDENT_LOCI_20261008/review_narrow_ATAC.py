from pathlib import Path
import sys,json,hashlib,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent
sys.path.insert(0,str(P/'06_v1.21.8/code'));sys.path.insert(1,str(P/'06_v1.21.1/code'))
import score_ship_candidates_v2 as legacy
import pyBigWig
cmd=json.loads((P/'06_v1.21.8/scoring_command.json').read_text());threshold=float(cmd[cmd.index('--min-atac-accessibility-percentile')+1]);track=Path(cmd[cmd.index('--atac-mean-bw')+1])
rows=json.loads((N/'next_queue.json').read_text());rows.append(dict(id='ANO2_39708',chrom='NC_051831.1',start0=39708074,end0=39709074))
bw=pyBigWig.open(str(track));bg=legacy.build_narrow_window_background(bw)
results=[]
for r in rows:
 s=int(r['start0']);e=int(r['end0']);mid=(s+e)//2;v=legacy.narrow_window_mean(bw,r['chrom'],s,e);p=legacy.percentile_rank(bg,v)
 results.append(dict(id=r['id'],chrom=r['chrom'],context_start0=max(0,mid-10000),context_end0=mid+10000,mean=v,percentile=p,passes_narrow_gate=p is not None and p>=threshold))
out=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),threshold=threshold,background_n=len(bg),background=bg.tolist(),track=str(track),track_sha256=hashlib.sha256(track.read_bytes()).hexdigest(),legacy_script_sha256=hashlib.sha256(Path(legacy.__file__).read_bytes()).hexdigest(),results=results,limitations=['This applies only the legacy 20kb center gate, not the full original SHIP-region wide gate.','Blood ATAC is not direct T-cell accessibility validation.'])
(N/'ATAC_narrow_gate_20261008.json').write_text(json.dumps(out,indent=2));print(json.dumps(results,indent=2))
