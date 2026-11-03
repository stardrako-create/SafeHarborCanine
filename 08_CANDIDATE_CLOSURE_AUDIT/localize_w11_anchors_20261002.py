from pathlib import Path
import pysam,json,datetime,hashlib
A=Path(__file__).resolve().parent; O=A/'W11_ANCHOR_LOCALIZATION_20261002';O.mkdir(exist_ok=False)
J=A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa'
with pysam.FastaFile(str(J)) as f:s=f.fetch('junction_L4707939_R33613308').upper()
assert len(s)==2000
patterns={**{'left_'+str(n):s[1000-n:1000] for n in [50,100,114]},**{'right_'+str(n):s[1000:1000+n] for n in [50,100,119,150,200,300]}}
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results={}
for label,ref in [('ROS',A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'),('UU',Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa'))]:
 rows={k:[] for k in patterns}
 with pysam.FastaFile(str(ref)) as f:
  for chrom in f.references:
   seq=f.fetch(chrom).upper()
   for k,p in patterns.items():
    for strand,q in [('+',p),('-',rc(p))]:
     start=0
     while True:
      i=seq.find(q,start)
      if i<0:break
      rows[k].append([chrom,i,i+len(q),strand]);start=i+1
 results[label]=rows
 (O/(label+'_exact_hits.json')).write_text(json.dumps(rows,indent=2))
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'patterns':patterns,'counts':{g:{k:len(v) for k,v in rows.items()} for g,rows in results.items()},'limits':['Exact matches only; mismatched copies may exist','Longer synthetic anchors beyond read coverage are not observed read evidence','Reference occurrences do not establish sample structural arrangement','Same sequence in distinct loci cannot localize source with these reads']}
(O/'review.json').write_text(json.dumps(summary,indent=2));(O/'script.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest());print(json.dumps(summary['counts']))