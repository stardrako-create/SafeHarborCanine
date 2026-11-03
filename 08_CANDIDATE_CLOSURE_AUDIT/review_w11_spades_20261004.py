from pathlib import Path
import pysam,json,hashlib,subprocess,datetime
A=Path(__file__).resolve().parent;O=A/'W11_SPADES_REVIEW_20261004';O.mkdir(exist_ok=False);C=Path('/mnt/d/SafeHarborWork/w11_spades_20261004_ascii/assembly_libfix/contigs.fasta')
assert hashlib.sha256(C.read_bytes()).hexdigest()==json.loads((A/'W11_SPADES_20261004/status_libfix.json').read_text())['sha256']
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';U=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with pysam.FastaFile(str(A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa')) as f:j=f.fetch('junction_L4707939_R33613308').upper()
with pysam.FastaFile(str(R)) as f:flank=rc(f.fetch('NC_051836.1',35342033,35342133).upper())
patterns={'first_junction_20_20':j[980:1020],'tail_junction_20_20':'A'*20+flank[:20]};rows=[]
with pysam.FastxFile(str(C)) as f:
 for x in f:
  seq=x.sequence.upper();hits={k:[{'strand':strand,'offset0':s.find(p)} for strand,s in [('+',seq),('-',rc(seq))] if p in s] for k,p in patterns.items()};rows.append({'name':x.name,'length':len(seq),'hits':hits})
(O/'motif_review.json').write_text(json.dumps({'contig_count':len(rows),'total_bases':sum(r['length'] for r in rows),'rows':rows,'limits':['Exact motif screen only; missing motif does not exclude variant sequence','Contigs require read support and genomic validation']},indent=2))
print(json.dumps({'contigs':len(rows),'hits':[r for r in rows if any(r['hits'].values())]}),flush=True)
for label,ref in [('ROS',R),('UU',U)]:
 (O/'status.json').write_text(json.dumps({'stage':'aligning','reference':label,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}))
 with (O/(label+'.sam')).open('w') as out,(O/(label+'.log')).open('w') as err:subprocess.run(['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','6',str(ref),str(C)],stdout=out,stderr=err,check=True)
(O/'status.json').write_text(json.dumps({'stage':'alignments_completed_needs_review','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}))