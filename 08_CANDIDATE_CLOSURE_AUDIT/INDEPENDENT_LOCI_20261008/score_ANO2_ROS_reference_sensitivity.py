from pathlib import Path
import json,pysam,hashlib,subprocess,math
N=Path(__file__).resolve().parent;P=N.parent.parent;S=N/'ANO2_conservation_identity_20261009_v3';O=N/'ANO2_ROS_model_sensitivity_20261009';O.mkdir(exist_ok=False)
m=json.loads((S/'base_mapping.json').read_text())['mapping'];lookup={v[1]:i for i,v in enumerate(m)};fa=pysam.FastaFile(str(P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'));chrom='NC_051831.1';start=39708074;ref=fa.fetch(chrom,start,start+1000).upper();seen={};parts=[];buf=[]
def consume(lines):
 rows=[l.split() for l in lines if l.startswith(('s ','s\t'))]
 if not rows:return
 dog=rows[0];assert dog[1]=='Canis_lupus_familiaris.chr27' and dog[4]=='+';pos=int(dog[2]);dna=[];offsets=[]
 for base in dog[6]:
  if base=='-':dna.append('-');continue
  if pos in lookup:
   offset=lookup[pos];offsets.append(offset);dna.append(ref[offset]);assert offset not in seen;seen[offset]=ref[offset]
  else:dna.append('-')
  pos+=1
 assert pos==int(dog[2])+int(dog[3])
 if not offsets:return
 assert offsets==list(range(min(offsets),max(offsets)+1))
 dog=[dog[0],'Canis_lupus_familiaris.'+chrom,str(start+min(offsets)),str(len(offsets)),'+',str(fa.get_reference_length(chrom)),''.join(dna)];rows[0]=dog;parts.append('a\n'+'\n'.join(' '.join(r) for r in rows)+'\n\n')
for line in (S/'ANO2_39708.maf').open():
 if not line.strip():consume(buf);buf=[]
 elif not line.startswith('#'):buf.append(line.rstrip())
consume(buf);assert set(seen)==set(range(1000)) and ''.join(seen[i] for i in range(1000))==ref
maf=O/'ROS_reference.maf';maf.write_text('##maf version=1\n\n'+''.join(parts));model=P/'06_v1.21.9/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod';wig=O/'ROS_reference.wig';cmd=['/home/stardrako/miniforge3/envs/cactus/bin/phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(maf)];(O/'command.json').write_text(json.dumps(cmd))
with wig.open('w') as out,(O/'phyloP.err').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
values={};pos=None
for line in wig.open():
 if line.startswith('fixedStep'):f=dict(x.split('=') for x in line.split()[1:]);pos=int(f['start'])-1;assert int(f.get('step',1))==int(f.get('span',1))==1
 elif line.strip() and not line.startswith(('#','track')):v=float(line);assert pos is not None and pos not in values and math.isfinite(v);values[pos]=v;pos+=1
assert set(values)==set(range(start,start+1000));scores=[values[i] for i in range(start,start+1000)];max50=max(sum(scores[i:i+50])/50 for i in range(951));original=json.loads((S/'receipt.json').read_text());result=dict(state='completed',reference_sequence_matches_ROS=True,scored_bp=1000,ROS_max50bp_mean=max50,CanFam3_transferred_max50bp_mean=original['max50bp_mean'],exceeds_legacy_6_5=max50>6.5,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'ANO2_39708.maf',S/'base_mapping.json',maf,wig,model]},method='Replace only dog reference with known ROS assembly sequence at 1000 unique monotone mapped columns; four base substitutions and 11 CanFam3-only bases represented as gaps in ROS; all other species alignment rows retained.',limitations=['Model sensitivity on retained alignment, not a new independent whole-genome alignment or biological validation.','No unsupported insertion columns invented; existing chain maps each ROS base uniquely.']);(O/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))