from pathlib import Path
import json,pysam,subprocess,hashlib
O=Path(__file__).parent
I=O/'allele_projection_20261006T115816Z/review.json'
D=O/'normalized_exact_blocks_20261006';D.mkdir(exist_ok=False)
R=O.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'
rows=[x for x in json.loads(I.read_text())['alleles'] if x['status']=='exact_block_REF_verified_normalization_pending']
lookup={};lines=[]
with pysam.FastaFile(str(R)) as fa:
 for i,x in enumerate(rows):
  chrom,start,end=x['window'].rsplit('_',2);pos=int(start)+x['ROS_offset0'];ident=f'a{i}'
  assert fa.fetch(chrom,pos,pos+len(x['ROS_ref'])).upper()==x['ROS_ref']
  lookup[ident]=dict(x,chrom=chrom,pos0=pos)
  lines.append((chrom,pos,f"{chrom}\t{pos+1}\t{ident}\t{x['ROS_ref']}\t{x['ROS_alt']}\t.\t.\t.\n"))
 header='##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c},length={fa.get_reference_length(c)}>\n' for c in sorted({x[0] for x in lines}))+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
 inp=D/'input.vcf';inp.write_text(header+''.join(x[2] for x in sorted(lines)))
 cmd=['/home/stardrako/miniforge3/envs/atac/bin/bcftools','norm','-f',str(R),'--check-ref','e','-Ov','-o',str(D/'normalized.vcf'),str(inp)]
 (D/'command.json').write_text(json.dumps(cmd));ret=subprocess.run(cmd,capture_output=True,text=True);(D/'norm.log').write_text(ret.stdout+ret.stderr);assert ret.returncode==0
 results=[]
 with pysam.VariantFile(str(D/'normalized.vcf')) as vf:
  for r in vf:
   old=lookup[r.id];lo=max(0,min(old['pos0'],r.start)-50);hi=max(old['pos0']+len(old['ROS_ref']),r.start+len(r.ref))+50
   seq=fa.fetch(r.contig,lo,hi).upper();a=seq[:old['pos0']-lo]+old['ROS_alt']+seq[old['pos0']-lo+len(old['ROS_ref']):];b=seq[:r.start-lo]+r.alts[0]+seq[r.start-lo+len(r.ref):]
   assert a==b
   results.append(dict(id=r.id,window=old['window'],chrom=r.contig,start0=r.start,ref=r.ref,alt=r.alts[0],source_allele=old['allele'],haplotype_equivalence_verified=True))
assert len(results)==len(rows)==2014
(D/'review.json').write_text(json.dumps(dict(results=results,input_sha256=hashlib.sha256(I.read_bytes()).hexdigest(),count=len(results),limitations=['Window-allele observations overlap; not 2014 unique population variants','422 observations remain unresolved','No candidate approval from normalization alone']),indent=2))
print('2014 allele observations normalized; REF and haplotype equivalence verified')