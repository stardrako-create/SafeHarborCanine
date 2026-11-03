from pathlib import Path
import pysam,json,datetime,hashlib
A=Path(__file__).resolve().parent;O=A/'W11_TAIL_REFERENCE_20261002';O.mkdir(exist_ok=False)
specs=[('PBGV000010',A/'SV_PBGV000010/alignment_chunked/selected_markdup.bam',Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa'),'NC_049253.1',[4707925,4707939],4707925,'left','A'),('PBGV000012',A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam',A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa','NC_051836.1',[35342119,35342133],35342133,'right','T')]
results={}
for sample,bam,ref,chrom,bounds,tailb,side,base in specs:
 counts={str(b):[] for b in bounds};tails=[]
 with pysam.AlignmentFile(str(bam),'rb') as f,pysam.FastaFile(str(ref)) as fa:
  for r in f.fetch(chrom,min(bounds)-200,max(bounds)+200):
   if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail or r.mapping_quality<30:continue
   mapping={p:q for q,p in r.get_aligned_pairs(matches_only=True)}
   for b in bounds:
    positions=list(range(b-20,b+20));qs=[mapping.get(p) for p in positions]
    if all(q is not None for q in qs) and qs==list(range(qs[0],qs[0]+40)) and min(r.query_qualities[q] for q in qs)>=20 and ''.join(r.query_sequence[q] for q in qs).upper()==fa.fetch(chrom,b-20,b+20).upper():counts[str(b)].append(r.query_name)
   op=r.cigartuples[0] if side=='left' else r.cigartuples[-1];boundary=r.reference_start if side=='left' else r.reference_end
   if op[0]!=4 or op[1]<20 or boundary!=tailb:continue
   n=op[1]
   if side=='left':clipqs=list(range(n-20,n));anchorpos=list(range(tailb,tailb+20))
   else:clipqs=list(range(len(r.query_sequence)-n,len(r.query_sequence)-n+20));anchorpos=list(range(tailb-20,tailb))
   aq=[mapping.get(p) for p in anchorpos];clipok=all(r.query_sequence[q]==base and r.query_qualities[q]>=20 for q in clipqs)
   anchorok=all(q is not None for q in aq) and min(r.query_qualities[q] for q in aq)>=20 and ''.join(r.query_sequence[q] for q in aq).upper()==fa.fetch(chrom,anchorpos[0],anchorpos[-1]+1).upper()
   tails.append({'read':r.query_name,'mate':1 if r.is_read1 else 2,'clip20_exact_BQ20':clipok,'reference_anchor20_exact_BQ20':anchorok,'passes_both':clipok and anchorok})
 results[sample]={'exact_reference_spanning_names':{b:sorted(set(v)) for b,v in counts.items()},'reference_counts':{b:len(set(v)) for b,v in counts.items()},'tail_rows':tails,'tail_pass_names':len({x['read'] for x in tails if x['passes_both']})}
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':results,'limits':['Exact20bp reference plus20bp homopolymer BQ20 is local-tail support, not full insertion or length','Reference spanning requires40 consecutive exact bases BQ20; counts not genotypes','Allelic bias, repeats and shared read evidence remain','No full allele connects the two boundaries yet']}
(O/'review.json').write_text(json.dumps(out,indent=2));(O/'script.sha256').write_text(hashlib.sha256(Path(__file__).read_bytes()).hexdigest());print(json.dumps({s:{'reference':v['reference_counts'],'tail_pass':v['tail_pass_names'],'tails':v['tail_rows']} for s,v in results.items()}))