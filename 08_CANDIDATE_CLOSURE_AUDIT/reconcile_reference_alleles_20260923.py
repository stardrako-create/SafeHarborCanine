from pathlib import Path
import json,pysam,collections,hashlib,datetime,subprocess
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';Q=C/'variant_catalogue_partial_query_20260923';O=C/'reference_allele_reconciliation_20260923';O.mkdir(exist_ok=False)
D=json.loads((C/'pending_allele_diagnostics_20260923/review.json').read_text());wanted={x['id']:x for x in D['results'] if 'ros_reference' in x};rows=[];rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
W={x['id']:x for x in json.loads((C/'gapped_projection_20260923/UU_review.json').read_text())['results']};R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';U='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa'
with pysam.FastaFile(str(R)) as ref,pysam.FastaFile(U) as uu:
 for file in sorted(Q.glob('bg*.vcf')):
  flat=0;window=file.name.split('__')[0];w=W[window]
  with pysam.VariantFile(str(file)) as vf:
   for record_index,r in enumerate(vf):
    ids=[file.stem+'__'+str(flat+i+1) for i in range(len(r.alts or ()))];flat+=len(ids);targets=[wanted[k] for k in ids if k in wanted]
    if not targets:continue
    d=targets[0];assert all(x['ros_start0']==d['ros_start0'] and x['ros_reference']==d['ros_reference'] for x in targets)
    original=[r.ref]+list(r.alts);alleles=[s.upper() if w['strand']=='+' else rc(s.upper()) for s in original];ros=d['ros_reference'];assert ref.fetch(d['ros_chrom'],d['ros_start0'],d['ros_end0']).upper()==ros
    result={'source_file':file.name,'source_record_index0':record_index,'pending_ids':[x['id'] for x in targets],'window':window,'ros_chrom':d['ros_chrom'],'ros_start0':d['ros_start0'],'ros_reference':ros,'source_alleles_oriented':alleles,'original_AF':list(r.info.get('AF',())),'filters':list(r.filter)}
    ac=r.info.get('AC');an=r.info.get('AN');ac=list(ac) if ac is not None else None
    result.update(original_AC=ac,original_AN=an)
    if ros not in alleles:result['status']='ROS_allele_absent_from_source_record';rows.append(result);continue
    if len(set(alleles))!=len(alleles) or any(not s or set(s)-set('ACGT') for s in alleles):result['status']='noncanonical_or_duplicate_alleles_need_context';rows.append(result);continue
    if ac is None or an is None or len(ac)!=len(alleles)-1 or any(x is None or x<0 for x in ac) or an<=0 or sum(ac)>an:result['status']='invalid_or_missing_AC_AN';rows.append(result);continue
    counts=[an-sum(ac)]+ac;ri=alleles.index(ros);order=[ri]+[i for i in range(len(alleles)) if i!=ri]
    assert sum(counts)==an and sorted(order)==list(range(len(alleles)))
    # Preserve complete original allele-count vector under a permutation, not an AF complement.
    result.update(status='allele_set_rebased_counts_preserved',source_index_for_ROS_REF=ri,source_indices_in_ROS_order=order,ROS_alleles=[alleles[i] for i in order],ROS_counts=[counts[i] for i in order],ROS_count_frequencies=[counts[i]/an for i in order],frequency_method='AC/AN; source REF count = AN - sum(all ALT AC); original rounded AF retained separately')
    assert sum(result['ROS_counts'])==an
    rows.append(result)
proposed={}
for i,r in enumerate(rows):
 if r['status']!='allele_set_rebased_counts_preserved':continue
 for j,alt in enumerate(r['ROS_alleles'][1:],1):proposed[f'R{i}_A{j}']={'record_index':i,'chrom':r['ros_chrom'],'start0':r['ros_start0'],'ref':r['ros_reference'],'alt':alt,'count':r['ROS_counts'][j],'AN':r['original_AN']}
with pysam.FastaFile(str(R)) as ref:
 header='##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c},length={ref.get_reference_length(c)}>\n' for c in sorted({r['chrom'] for r in proposed.values()}))+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
 (O/'rebased_before_norm.vcf').write_text(header+''.join(f"{r['chrom']}\t{r['start0']+1}\t{k}\t{r['ref']}\t{r['alt']}\t.\t.\t.\n" for k,r in proposed.items()))
 cmd=['/home/stardrako/miniforge3/envs/atac/bin/bcftools','norm','-f',str(R),'--check-ref','e','-m','-any','-Ov','-o',str(O/'rebased_normalized.vcf'),str(O/'rebased_before_norm.vcf')]
 p=subprocess.run(cmd,capture_output=True,text=True);(O/'norm.log').write_text(p.stderr);assert p.returncode==0
 normalized=[];seen=set()
 with pysam.VariantFile(str(O/'rebased_normalized.vcf')) as vf:
  for r in vf:
   assert r.id not in seen;seen.add(r.id);v=proposed[r.id];p=v['start0'];assert ref.fetch(r.contig,r.start,r.stop).upper()==r.ref.upper()
   lo=max(0,min(p,r.start)-50);hi=min(ref.get_reference_length(r.contig),max(p+len(v['ref']),r.stop)+50)
   assert ref.fetch(r.contig,lo,p).upper()+v['alt']+ref.fetch(r.contig,p+len(v['ref']),hi).upper()==ref.fetch(r.contig,lo,r.start).upper()+r.alts[0].upper()+ref.fetch(r.contig,r.stop,hi).upper()
   normalized.append(dict(v,id=r.id,start0=r.start,ref=r.ref.upper(),alt=r.alts[0].upper()))
 assert seen==set(proposed)
assert {k for r in rows for k in r['pending_ids']}==set(wanted)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pending_allele_entries_examined':len(wanted),'source_records_examined':len(rows),'record_status_counts':dict(collections.Counter(r['status'] for r in rows)),'normalized_rebased_ALT_records':len(normalized),'records':rows,'limitations':['Locus allele-set representation only, not whole-window clearance or donor genotyping','Counts are permuted from each source record, never pooled across records or treated as phase','Normalized ROS haplotype equivalence checked; cross-assembly local sequence differences outside the locus remain','Kept separate from 319 full-window and 334 partial exact-run allele records','Not a final solution for symbolic alleles, mapping gaps or w11 structural uncertainty']}
(O/'review.json').write_text(json.dumps(summary,indent=2));(O/'normalized_alleles.json').write_text(json.dumps(normalized,indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ['records','limitations']},indent=2))
