from pathlib import Path
import json,pysam,collections,hashlib
O=Path(__file__).parent
maps={x['id']:x for x in json.loads((O/'basewise_projection.json').read_text())['results']}
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results=[]
with pysam.FastaFile(str(O/'UU_projection/controls_ROS.fa')) as ros:
 for folder in ['variant_catalogue_query_drop_samples','variant_catalogue_query_remaining3']:
  for file in sorted((O/folder).glob('*__*.vcf')):
   w=file.name.split('__')[0];m=maps[w];seq=ros.fetch(w).upper();mm={v['UU_position0']:v for v in m['mismatches']}
   with pysam.VariantFile(str(file)) as vf:
    for r in vf:
     if len(r.ref)!=1 or r.start not in mm:continue
     alleles=[r.ref]+list(r.alts or ())
     if any(len(a)!=1 or a.upper() not in 'ACGT' for a in alleles):continue
     oriented=[a.upper() if m['strand']=='+' else rc(a.upper()) for a in alleles]
     offset=mm[r.start]['ROS_offset0'];ref=seq[offset]
     if oriented.count(ref)!=1:continue
     ac=list(r.info.get('AC',()));an=r.info.get('AN');assert len(ac)==len(alleles)-1 and an and 0<=sum(ac)<=an
     counts=[an-sum(ac)]+ac;index=oriented.index(ref);order=[index]+[i for i in range(len(alleles)) if i!=index]
     reordered=[counts[i] for i in order];assert sum(reordered)==an
     results.append(dict(window=w,ROS_offset0=offset,ROS_alleles=[oriented[i] for i in order],ROS_counts=reordered,AN=an,ROS_frequencies=[counts[i]/an for i in order],source_REF=r.ref,source_ALT=list(r.alts),source_AC=ac,source_AF=list(r.info.get('AF',())),source_file=str(file),source_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),filters=list(r.filter),source_indices=order))
report=dict(records=results,count=len(results),method='Permutation of full allele set; source REF count AN minus sum of all AC. ROS REF checked against extracted reference.',limitations=['No population counts pooled across overlapping windows','Not individual genotype or haplotype validation','Only single-base sequence allele sets covered; symbolic/gapped cases remain pending'])
(O/'SNP_reference_rebase_counts_20261008.json').write_text(json.dumps(report,indent=2))
print(json.dumps({'records':len(results),'priority':[r for r in results if r['window']=='NC_051807.1_2091307_2092307']}))
