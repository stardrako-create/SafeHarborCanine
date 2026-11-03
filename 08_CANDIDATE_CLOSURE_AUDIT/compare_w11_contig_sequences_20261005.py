from pathlib import Path
import pysam,json
A=Path(__file__).resolve().parent;out={};seqs=[]
for name,path,node in [('PBGV10','/mnt/d/SafeHarborWork/w11_PBGV10_20261005/assembly/contigs.fasta','NODE_22_length_154_cov_6.896104'),('PBGV12','/mnt/d/SafeHarborWork/w11_spades_20261004_ascii/assembly_libfix/contigs.fasta','NODE_20_length_175_cov_6.704082')]:
 with pysam.FastxFile(path) as f:
  seq=next(x.sequence for x in f if x.name==node)
 seq=seq.translate(str.maketrans('ACGT','TGCA'))[::-1];seqs.append(seq);out[name]={'node':node,'oriented_sequence':seq,'junction_offset0':77}
common=0
for a,b in zip(*seqs):
 if a!=b:break
 common+=1
out['common_exact_prefix_bp']=common;out['first_is_prefix_of_second']=seqs[1].startswith(seqs[0]);out['limits']=['Same assembly workflow on distinct samples; not orthogonal validation','Shared partial sequence does not resolve full allele or distant origin'];(A/'W11_SPADES_PBGV10_REVIEW_20261005/cross_sample_junction.json').write_text(json.dumps(out,indent=2));print(common,out['first_is_prefix_of_second'])