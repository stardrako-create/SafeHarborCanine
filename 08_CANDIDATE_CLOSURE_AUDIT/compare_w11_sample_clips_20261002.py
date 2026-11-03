from pathlib import Path
import pysam,json,datetime,collections,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_TWO_SAMPLE_CLIPS_20261002';O.mkdir(exist_ok=False)
B=A/'SV_PBGV000010/alignment_chunked/selected_markdup.bam';rows=[]
with pysam.AlignmentFile(str(B),'rb') as f:
 chrom='NC_049253.1'
 assert chrom in f.references
 for r in f.fetch(chrom,4707400,4708450):
  if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality<30:continue
  for side,op in [('left',r.cigartuples[0]),('right',r.cigartuples[-1])]:
   if op[0]!=4 or op[1]<20:continue
   n=op[1];seq=r.query_sequence[:n] if side=='left' else r.query_sequence[-n:];q=list(r.query_qualities[:n] if side=='left' else r.query_qualities[-n:])
   rows.append({'name':r.query_name,'mate':1 if r.is_read1 else 2,'side':side,'boundary0':r.reference_start if side=='left' else r.reference_end,'clip':seq,'all_BQ20':min(q)>=20,'BQ20_bases':sum(x>=20 for x in q),'MAPQ':r.mapping_quality,'cigar':r.cigarstring})
groups=[]
for side,b in sorted({(x['side'],x['boundary0']) for x in rows}):
 rr=[x for x in rows if x['side']==side and x['boundary0']==b];groups.append({'side':side,'boundary0':b,'reads':len(rr),'all_BQ20':sum(x['all_BQ20'] for x in rr),'clips':[x['clip'] for x in rr]})
ref=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');ros=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with pysam.FastaFile(str(ref)) as f,pysam.FastaFile(str(ros)) as g:
 target=rc(g.fetch('NC_051836.1',35342119,35342133).upper());local=f.fetch('NC_049253.1',4707400,4708450).upper();positions=[];start=0
 while True:
  p=local.find(target,start)
  if p<0:break
  positions.append(4707400+p);start=p+1
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sample':'PBGV000010','reference':'UU','groups':groups,'rows':rows,'ROS14bp_reverse_sequence':target,'exact14bp_occurrences_in_UU_window0':positions,'limits':['Same short-read evidence, not orthogonal validation','Homopolymer clips and14bp spacing do not prove insertion length or full allele','Reference-local14bp matching is not rigorous whole-flank projection']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'script.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest());print(json.dumps({'groups':groups,'positions':positions}))