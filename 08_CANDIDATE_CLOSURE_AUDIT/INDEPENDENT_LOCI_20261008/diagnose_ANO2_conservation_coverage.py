from pathlib import Path
import json,math,hashlib,collections
N=Path(__file__).resolve().parent;P=N.parent.parent;mapping=json.loads((N/'ANO2_conservation_identity_20261009_v2/base_mapping.json').read_text());m=mapping['mapping'];want={v[1] for v in m};source=Path('/mnt/d/Jin2024_work/zoonomia_hal/maf_dog/cand_36.chr27.39320764.maf');covered=set();duplicates=set()
for line in source.open():
 f=line.split()
 if len(f)==7 and f[0]=='s' and f[1]=='Canis_lupus_familiaris.chr27':
  assert f[4]=='+';s=int(f[2]);n=int(f[3]);assert len(f[6].replace('-',''))==n
  q=set(range(s,s+n))&want;duplicates.update(q&covered);covered.update(q)
wig=P/'06_v1.21.9/ANO2.published_model.wig';values={};pos=None
for line in wig.open():
 if line.startswith('fixedStep'):f=dict(x.split('=') for x in line.split()[1:]);pos=int(f['start'])-1;step=int(f.get('step',1))
 elif line.strip() and not line.startswith(('#','track')):values[pos]=float(line);pos+=step
missing=[dict(ROS_offset0=i,ROS_position0=39708074+i,CanFam3_position0=v[1],in_MAF=v[1] in covered,in_WIG=v[1] in values,WIG_value=values.get(v[1])) for i,v in enumerate(m) if v[1] not in covered or v[1] not in values or not math.isfinite(values[v[1]])]
result=dict(reference_mapped_bp=1000,MAF_reference_covered_bp=len(covered),WIG_finite_mapped_bp=sum(v[1] in values and math.isfinite(values[v[1]]) for v in m),duplicate_reference_positions=sorted(duplicates),missing=missing,chain_jumps=mapping['jumps'],envelope_bp=mapping['envelope_bp'],historical_claim_scored_fraction=1.0,interpretation='Historic scored_fraction is not accepted until reconciled against basewise source coordinates. No missing position receives zero or PASS.',source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,wig]});out=N/'ANO2_conservation_coverage_discrepancy_20261009.json';assert not out.exists();out.write_text(json.dumps(result,indent=2,allow_nan=False));print(json.dumps(result,allow_nan=False))