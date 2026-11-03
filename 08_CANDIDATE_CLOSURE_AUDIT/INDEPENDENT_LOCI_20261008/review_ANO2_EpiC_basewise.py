from pathlib import Path
import json,gzip,bisect,hashlib,collections,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent;m=json.loads((N/'ANO2_conservation_identity_20261009_v3/base_mapping.json').read_text())['mapping'];positions=sorted(v[1] for v in m);assert len(set(positions))==1000;chrom=m[0][0];files=sorted((P/'06_v1.21.5/external_sources').glob('*_13_dense.bed.gz'));assert files;results=[]
for file in files:
 counts=collections.Counter();labels={};hits=[]
 with gzip.open(file,'rt') as f:
  for line in f:
   if line.startswith(('#','track','browser')):continue
   x=line.split()
   if len(x)<4 or x[0]!=chrom:continue
   s=int(x[1]);e=int(x[2]);state=int(x[3]);lo=bisect.bisect_left(positions,s);hi=bisect.bisect_left(positions,e)
   if hi<=lo:continue
   for pos in positions[lo:hi]:assert pos not in labels;labels[pos]=state
   counts[state]+=hi-lo
   if 1<=state<=7:hits.append(dict(start0=s,end0=e,state=state,ROS_target_bp=hi-lo))
 results.append(dict(source=str(file),tissue=file.name.split('_')[0],sha256=hashlib.sha256(file.read_bytes()).hexdigest(),mapped_ROS_bp=1000,annotated_ROS_bp=len(labels),state_bp=dict(counts),promoter_enhancer_ROS_bp=sum(v for k,v in counts.items() if 1<=k<=7),promoter_enhancer_hits=hits))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),window='NC_051831.1_39708074_39709074',mapping_sha256=hashlib.sha256((N/'ANO2_conservation_identity_20261009_v3/base_mapping.json').read_bytes()).hexdigest(),results=results,complete_annotation_in_all_sources=all(r['annotated_ROS_bp']==1000 for r in results),any_promoter_enhancer_hit=any(r['promoter_enhancer_ROS_bp'] for r in results),limitations=['Exact per-base intersection of all1000 uniquely projected ROS positions, not a contiguous envelope approximation.','Source state1-7 interpretation unchanged from existing pipeline.','No CAR-T-specific regulatory isolation or biological safety demonstration.']);out=N/'ANO2_EpiC_basewise_review_20261009.json';assert not out.exists();out.write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='results'}));print(json.dumps([dict(tissue=r['tissue'],coverage=r['annotated_ROS_bp'],PE=r['promoter_enhancer_ROS_bp']) for r in results]))