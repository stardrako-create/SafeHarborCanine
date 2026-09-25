from pathlib import Path
import json,csv,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.2';OLD=P/'06_v1.22.1'
m=next(r for r in csv.DictReader((OLD/'dog10k_mapping.tsv').open(),delimiter='\t') if r['window_id']=='w00');s=int(m['start']);e=int(m['end']);rs=int(m['ros_start']);fa=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa')
with pysam.AlignmentFile(str(OLD/'windows_uu.sam'),'r') as f:a=next(a for a in f if a.query_name=='w00' and not a.is_secondary and not a.is_supplementary)
assert not a.is_reverse;project={r:rs+q for q,r in a.get_aligned_pairs() if q is not None and r is not None};out=[];summ=[];base='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/'
for name in ['AutoAndXPAR.SNPs.vqsr99.vcf.gz','AutoAndXPAR.nonSNPs.filter.vcf.gz']:
 vf=pysam.VariantFile(base+name);rows=list(vf.fetch(m['vcf_contig'],s,e));rows=[r for r in rows if r.start<e and r.stop>s]
 with (O/(name+'.LOC_supplementary.vcf')).open('w') as f:f.write(str(vf.header));f.writelines(str(r) for r in rows)
 common=0;passed=0;uncertain=0
 for r in rows:
  assert fa.fetch(m['uu_contig'],r.start,r.start+len(r.ref)).upper()==r.ref.upper()
  af=r.info.get('AF');afs=af if isinstance(af,tuple) else (af,);mx=max((x for x in afs if x is not None),default=None);isp=list(r.filter)==['PASS'];passed+=isp;common+=isp and mx is not None and mx>=.01
  positions=[project.get(p) for p in range(r.start,r.start+len(r.ref))];exact=all(p is not None for p in positions) and all(b==a+1 for a,b in zip(positions,positions[1:]));uncertain+=not exact
  out.append(dict(callset=name,uu_chrom=r.chrom,uu_start=r.start,uu_end=r.stop,REF=r.ref,ALT=','.join(r.alts or []),FILTER=';'.join(r.filter),AF=json.dumps(af),max_alt_AF=mx,AC=json.dumps(r.info.get('AC')),AN=r.info.get('AN'),REF_matches_UU=True,REF_span_projects_contiguously=exact,ros_REF_span_start=positions[0] if exact else '',ros_REF_span_end=positions[-1]+1 if exact else '',projection_scope='REF-span coordinates only; ALT not normalized onto ROS'))
 summ.append(dict(callset=name,samples=len(vf.header.samples),records=len(rows),PASS=passed,PASS_AF_ge_0_01=common,noncontiguous_or_missing_REF_projection=uncertain,url=base+name));vf.close();print(name,len(rows),flush=True)
with (O/'LOC_supplementary_variants.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
sv=pysam.VariantFile('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz');svs=[dict(id=v.id,start=v.start,end=v.stop,filter=';'.join(v.filter)) for v in sv.fetch(m['vcf_contig'],s,e) if v.start<e and v.stop>s]
(O/'LOC_supplementary_lookup.json').write_text(json.dumps(dict(window_id='w00',UU_span=[m['vcf_contig'],s,e],summary=summ,SV_records=svs,limitations=['Supplementary UU reference-span lookup; original automated NM gate still failed','Interpret alongside reciprocal mapping audit','16 ROS-specific unpaired bases cannot be assessed for SNPs from this UU VCF','Even contiguous REF-coordinate projection does not validate transferred ALT allele representation','No donor genotype inference or safety nomination']),indent=2));print(json.dumps(summ,indent=2))
