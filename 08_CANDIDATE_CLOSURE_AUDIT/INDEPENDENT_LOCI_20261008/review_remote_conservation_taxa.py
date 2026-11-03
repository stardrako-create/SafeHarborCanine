from pathlib import Path
import json,collections,subprocess,math,hashlib,statistics
N=Path(__file__).resolve().parent;S=N/'remote_conservation_verified_20261009_v3';st=json.loads((S/'status.json').read_text());assert st['state']=='completed' and len(st['results'])==3 and not st['errors'];O=N/'remote_conservation_taxa_review_20261009';O.mkdir(exist_ok=False);model=N.parent.parent/'06_v1.21.9/Anc239_allARs_100kb_lessGC40_241species_30Consensus.mod';results=[]
for record in st['results']:
 name=record['id'];source=S/(name+'.maf');blocks=[];buf=[]
 for line in source.read_text().splitlines():
  if not line.strip():
   if buf:blocks.append(buf);buf=[]
  elif not line.startswith('#'):buf.append(line)
 if buf:blocks.append(buf)
 output=[];positions=[];counts=[];details=[];removed=collections.Counter()
 for block in blocks:
  rows=[l.split() for l in block if l.startswith(('s ','s\t'))];assert len([r for r in rows if r[1].startswith('Canis_lupus_familiaris.')])==1 and rows[0][1].startswith('Canis_lupus_familiaris.');assert rows[0][4]=='+';assert len({len(r[6]) for r in rows})==1;groups=collections.Counter(r[1].split('.',1)[0] for r in rows);dup={g for g,c in groups.items() if c>1};assert 'Canis_lupus_familiaris' not in dup
  keep=[r for r in rows if r[1].split('.',1)[0] not in dup];assert len({r[1].split('.',1)[0] for r in keep})==len(keep)
  for r in rows:assert sum(b!='-' for b in r[6])==int(r[3]) and 0<=int(r[2]) and int(r[2])+int(r[3])<=int(r[5]);assert set(r[6].upper())<=set('ACGTRYSWKMBDHVN-')
  a=[l for l in block if l.split() and l.split()[0]=='a'];assert len(a)==1;output.append('\n'.join(a+['\t'.join(r) for r in keep])+'\n\n');p=int(rows[0][2])
  for col,base in enumerate(rows[0][6]):
   if base=='-':continue
   assert base.upper() in 'ACGT';positions.append(p);p+=1;counts.append(sum(r[6][col].upper() in 'ACGT' for r in keep[1:]))
  for g in dup:removed[g]+=int(rows[0][3])
  details.append(dict(start0=int(rows[0][2]),bp=int(rows[0][3]),raw_rows=len(rows),unique_retained_rows=len(keep),removed_species=sorted(dup)))
 assert len(positions)==1000 and len(set(positions))==1000;assert sorted(positions)==list(range(min(positions),min(positions)+1000))
 maf=O/(name+'_unique_taxa.maf');maf.write_text('##maf version=1\n\n'+''.join(output));wig=O/(name+'_unique_taxa.wig');cmd=['/home/stardrako/miniforge3/envs/cactus/bin/phyloP','--method','LRT','--mode','CONACC','--wig-scores',str(model),str(maf)];(O/(name+'.command.json')).write_text(json.dumps(cmd))
 with wig.open('w') as out,(O/(name+'.err')).open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=600)
 vals={};p=None
 for line in wig.read_text().splitlines():
  if line.startswith('fixedStep'):
   v=dict(x.split('=') for x in line.split()[1:]);p=int(v['start'])-1;assert v.get('step','1')=='1' and v.get('span','1')=='1'
  elif line.strip() and not line.startswith(('#','track')):
   v=float(line);assert p is not None and p not in vals and math.isfinite(v);vals[p]=v;p+=1
 assert set(vals)==set(positions);vv=[vals[p] for p in sorted(vals)];max50=max(sum(vv[i:i+50])/50 for i in range(951));results.append(dict(id=name,state='completed',scored_reference_bp=1000,max50bp_mean=max50,exceeds_legacy_6_5=max50>6.5,raw_max50bp_mean=record['max50bp_mean'],removed_duplicate_species_refbp=dict(removed),nonreference_ACGT=dict(min=min(counts),median=statistics.median(counts),max=max(counts),reference_bp_with_zero_nonreference_ACGT=sum(c==0 for c in counts)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),curated_maf_sha256=hashlib.sha256(maf.read_bytes()).hexdigest(),wig_sha256=hashlib.sha256(wig.read_bytes()).hexdigest(),block_details=details))
receipt=dict(state='completed',results=results,source_receipt_sha256=hashlib.sha256((S/'status.json').read_bytes()).hexdigest(),policy='Remove all rows of any nonreference species duplicated within a block; retain original coordinates/DOGreference and unchanged phyloPmodel6.5gate; not scored-based selection.',limits=['Same HALalignment, not independent biologicalconfirmation.','Missing/ambiguous species data are not absence of conservation; informationcounts explicit.','This receipt assesses conservation only, not other locus/guide/experimentalcriteria.']);(O/'receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps([{k:v for k,v in r.items() if k!='block_details'} for r in results]))
