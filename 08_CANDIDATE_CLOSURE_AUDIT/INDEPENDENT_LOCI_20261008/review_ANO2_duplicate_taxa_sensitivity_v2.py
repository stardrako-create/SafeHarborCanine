from pathlib import Path
import json,collections,subprocess,hashlib,math,statistics
N=Path(__file__).resolve().parent;S=N/'ANO2_ROS_model_sensitivity_20261009';O=N/'ANO2_duplicate_taxa_sensitivity_20261009_v2';O.mkdir(exist_ok=False);M=N.parent.parent/'06_v1.21.9/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod';source=S/'ROS_reference.maf';blocks=[];buf=[]
for line in source.read_text().splitlines():
 if not line.strip():
  if buf:blocks.append(buf);buf=[]
 elif not line.startswith('#'):buf.append(line)
if buf:blocks.append(buf)
outputs={'duplicate_taxa_masked':[],'duplicate_rows_reversed':[]};details=[];counts=[];positions=[];consistency=collections.Counter()
for block in blocks:
 rows=[l.split() for l in block if l.startswith(('s ','s\t'))];assert len([r for r in rows if r[1].startswith('Canis_lupus_familiaris.')])==1;assert rows[0][1].startswith('Canis_lupus_familiaris.');group=collections.defaultdict(list)
 for row in rows:group[row[1].split('.',1)[0]].append(row)
 duplicates={k:v for k,v in group.items() if len(v)>1};assert 'Canis_lupus_familiaris' not in duplicates
 masked=[r for r in rows if r[1].split('.',1)[0] not in duplicates];reversed_rows=[];used=set()
 for r in rows:
  k=r[1].split('.',1)[0]
  if k in used:continue
  used.add(k);reversed_rows.extend(reversed(group[k]))
 prefix=[l for l in block if l.split() and l.split()[0]=='a'];assert len(prefix)==1
 for label,rr in [('duplicate_taxa_masked',masked),('duplicate_rows_reversed',reversed_rows)]:outputs[label].append('\n'.join(prefix+['\t'.join(r) for r in rr])+'\n\n')
 ref=rows[0];p=int(ref[2]);size=len(ref[6]);assert len({len(r[6]) for r in rows})==1
 for k,rr in duplicates.items():
  for col,b in enumerate(ref[6]):
   if b=='-':continue
   bases={r[6][col].upper() for r in rr if r[6][col].upper() in 'ACGT'}
   consistency['duplicate_species_refcolumns']+=1;consistency['duplicate_conflicting_ACGT_refcolumns']+=len(bases)>1
 for col,b in enumerate(ref[6]):
  if b=='-':continue
  assert b.upper() in 'ACGT';positions.append(p);p+=1;counts.append(sum(r[6][col].upper() in 'ACGT' for r in masked[1:]))
 details.append(dict(ref_start0=int(ref[2]),ref_bp=int(ref[3]),raw_rows=len(rows),retained_rows=len(masked),removed_species={k:[r[1] for r in rr] for k,rr in duplicates.items()}))
assert len(counts)==1000 and len(set(positions))==1000 and set(positions)==set(range(39708074,39709074))
def wig_values(path):
 values={};pos=None
 for line in path.read_text().splitlines():
  if line.startswith('fixedStep'):
   fields=dict(t.split('=') for t in line.split()[1:]);pos=int(fields['start'])-1;assert fields.get('step','1')=='1' and fields.get('span','1')=='1'
  elif line.strip() and not line.startswith(('#','track')):
   v=float(line);assert math.isfinite(v) and pos is not None and pos not in values;values[pos]=v;pos+=1
 assert set(values)==set(positions);return [values[p] for p in sorted(values)]
original=wig_values(S/'ROS_reference.wig');results=[]
for label,parts in outputs.items():
 maf=O/(label+'.maf');maf.write_text('##maf version=1\n\n'+''.join(parts));wig=O/(label+'.wig');cmd=['/home/stardrako/miniforge3/envs/cactus/bin/phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(M),str(maf)];(O/(label+'.command.json')).write_text(json.dumps(cmd))
 with wig.open('w') as out,(O/(label+'.err')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
 vals=wig_values(wig);max50=max(sum(vals[i:i+50])/50 for i in range(951));results.append(dict(label=label,scored_ROS_bp=1000,max50bp_mean=max50,exceeds_legacy_6_5=max50>6.5,changed_scores_vs_original=sum(a!=b for a,b in zip(vals,original)),max_absolute_change_vs_original=max(abs(a-b) for a,b in zip(vals,original)),maf_sha256=hashlib.sha256(maf.read_bytes()).hexdigest(),wig_sha256=hashlib.sha256(wig.read_bytes()).hexdigest()))
result=dict(state='completed',policy='All rows of any nonreference species with multiple rows in a block removed inmaskedanalysis; not selected by score or preferredsequence. Second analysis reverses duplicate roworder to characterize parser sensitivity.',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),model_sha256=hashlib.sha256(M.read_bytes()).hexdigest(),results=results,duplicate_consistency=dict(consistency),masked_nonreference_ACGT=dict(min=min(counts),median=statistics.median(counts),max=max(counts),zero_refbases=sum(c==0 for c in counts)),block_details=details,limits=['Sensitivity of same retainedalignment, not independent biological/evolutionarydata.','Removing ambiguous species rows reduces information; missingtaxa are not unconservedDNA.','No new speciescount or conservationthreshold imposed; original6.5gate unchanged.']);(O/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='block_details'}))

