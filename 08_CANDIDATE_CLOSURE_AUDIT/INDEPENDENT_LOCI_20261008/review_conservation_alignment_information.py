from pathlib import Path
import json,collections,hashlib,statistics,sys
N=Path(__file__).resolve().parent;include_remote='--remote' in sys.argv;datasets=[('ANO2_ROS',N/'ANO2_ROS_model_sensitivity_20261009/ROS_reference.maf',1000)]
if include_remote:
 S=N/'remote_conservation_verified_20261009_v3';st=json.loads((S/'status.json').read_text());assert st['state']=='completed' and not st['errors'];datasets += [(r['id'],S/(r['id']+'.maf'),1000) for r in st['results']]
results=[]
for name,path,expected in datasets:
 counts=[];ambiguous_nonreference=[];ref_positions=[];taxa=set();blocks=[];buf=[];duplicates=[];bad=[]
 def consume(lines):
  rows=[l.split() for l in lines if l.startswith(('s ','s\t'))]
  if not rows:return
  assert len({len(r[6]) for r in rows})==1
  refs=[r for r in rows if r[1].startswith('Canis_lupus_familiaris.')];assert len(refs)==1;ref=refs[0];assert ref[4]=='+';species=[r[1].split('.',1)[0] for r in rows];dup=[s for s,n in collections.Counter(species).items() if n>1];duplicates.extend(dup)
  for r in rows:
   assert sum(b!='-' for b in r[6])==int(r[3]);assert 0<=int(r[2]) and int(r[2])+int(r[3])<=int(r[5]);assert set(r[6].upper())<=set('ACGTRYSWKMBDHVN-');taxa.add(r[1].split('.',1)[0])
  others=[r for r in rows if r is not ref];pos=int(ref[2]);block_scores=[]
  for col,base in enumerate(ref[6]):
   if base=='-':continue
   assert base.upper() in 'ACGT';informative=sum(r[6][col].upper() in 'ACGT' for r in others);ambiguities=sum(r[6][col].upper() in 'RYSWKMBDHVN' for r in others);counts.append(informative);ambiguous_nonreference.append(ambiguities);ref_positions.append(pos);pos+=1;block_scores.append(informative)
  blocks.append(dict(reference_src=ref[1],start0=int(ref[2]),ref_bp=int(ref[3]),rows=len(rows),nonreference_taxa=len(others),min_nonreference_ACGT_per_ref_base=min(block_scores) if block_scores else None))
 for line in path.read_text().splitlines():
  if not line.strip():consume(buf);buf=[]
  elif not line.startswith('#'):buf.append(line)
 consume(buf);assert len(counts)==expected and len(set(ref_positions))==expected;assert not duplicates,('Duplicated species rows',duplicates)
 results.append(dict(id=name,source=str(path),source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),reference_bp=len(counts),blocks=len(blocks),taxa=len(taxa),min_nonreference_ACGT=min(counts),median_nonreference_ACGT=statistics.median(counts),max_nonreference_ACGT=max(counts),reference_bp_with_zero_nonreference_ACGT=sum(v==0 for v in counts),nonreference_ACGT_counts_per_base=counts,ambiguous_nonreference_counts_per_base=ambiguous_nonreference,block_details=blocks))
out=N/('conservation_alignment_information_review_20261009_v2.json' if include_remote else 'conservation_alignment_information_review_20261009.json');assert not out.exists();out.write_text(json.dumps(dict(state='completed',results=results,limits=['Alignment-information diagnostic, no new species-count acceptance threshold imposed.','Ambiguous/missing nonreference bases are missing information, not unconserved DNA.','Same underlyingalignment, not independent biologicalconfirmation.','WIG chrom label is inputfile stem; reference genomicidentity/coordinates checked viaMAF and independentDNA.']),indent=2));print(json.dumps([{k:v for k,v in r.items() if k not in ['nonreference_ACGT_counts_per_base','ambiguous_nonreference_counts_per_base','block_details']} for r in results]))
