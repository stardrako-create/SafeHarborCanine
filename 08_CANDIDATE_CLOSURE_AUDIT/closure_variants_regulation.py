from pathlib import Path
import csv,json,pysam,hashlib,subprocess,urllib.parse
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'08_CANDIDATE_CLOSURE_AUDIT';O.mkdir(exist_ok=True)
R=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';fa=pysam.FastaFile(str(R));uu=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');rows=list(csv.DictReader((P/'07_FINAL_CANDIDATES_2026-09-19/shortlist_evidence.tsv').open(),delimiter='\t'));variants=list(csv.DictReader((P/'06_v1.22.1/small_variant_records.tsv').open(),delimiter='\t'))
rc=lambda s:s.translate(str.maketrans('ACGTNacgtn','TGCANtgcan'))[::-1]
def write(name,rows):
 with (O/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
projected=[];vcf=[]
for w in rows:
 s,e=int(w['start']),int(w['end']);us,ue=int(w['uu_start']),int(w['uu_end']);strand=w['uu_strand'];assert w['uu_cigar']=='1000M' and w['uu_NM']=='0'
 a=fa.fetch(w['chrom'],s,e).upper();b=uu.fetch(w['uu_uu_contig'],us,ue).upper();assert a==(b if strand=='+' else rc(b))
 for i,v in enumerate(x for x in variants if x['window_id']==w['window_id']):
  vs,ve=int(v['start']),int(v['end']);ref=v['ref'].upper();alts=v['alt'].upper().split(',');assert uu.fetch(w['uu_uu_contig'],vs,vs+len(ref)).upper()==ref
  if not us<=vs<ve<=ue:raise RuntimeError('Variant extends beyond verified exact mapping')
  pos=s+vs-us if strand=='+' else e-(vs-us+len(ref));rr=ref if strand=='+' else rc(ref);aa=alts if strand=='+' else [rc(x) for x in alts]
  assert all(set(x)<=set('ACGTN') for x in [rr]+aa);assert fa.fetch(w['chrom'],pos,pos+len(rr)).upper()==rr
  ident=w['window_id']+'_'+str(i);projected.append(dict(id=ident,window_id=w['window_id'],UU_start0=vs,UU_ref=ref,UU_alt=v['alt'],ROS_chrom=w['chrom'],ROS_start0=pos,ROS_ref=rr,ROS_alt=','.join(aa),filter=v['filter'],AF=v['AF'],max_alt_AF=v['max_alt_AF'],projection='exact_full_window_sequence_verified; normalization_in_separate_VCF'))
  vcf.append((w['chrom'],pos,ident,rr,','.join(aa),v['filter']))
write('variant_projection_provenance.tsv',projected)
header='##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c},length={fa.get_reference_length(c)}>\n' for c in sorted({v[0] for v in vcf}))+''.join(f'##FILTER=<ID={x},Description="Original catalogue FILTER">\n' for x in sorted({f for v in vcf for f in v[5].split(';')} - {'PASS','.',''}))+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
inp=O/'projected_unnormalized.vcf';inp.write_text(header+''.join(f'{c}\t{pos+1}\t{i}\t{r}\t{a}\t.\t{f}\t.\n' for c,pos,i,r,a,f in sorted(vcf)))
cmd=['/home/stardrako/miniforge3/envs/atac/bin/bcftools','norm','-f',str(R),'--check-ref','e','-m','-any','-Ov','-o',str(O/'projected_normalized_ROS.vcf'),str(inp)]
ret=subprocess.run(cmd,capture_output=True,text=True);(O/'bcftools_norm.log').write_text(ret.stdout+ret.stderr);assert ret.returncode==0
lookup={r['id']:r for r in projected};norm=[]
for r in pysam.VariantFile(str(O/'projected_normalized_ROS.vcf')):
 rr=r.ref.upper();alt=r.alts[0].upper()
 assert fa.fetch(r.chrom,r.start,r.start+len(rr)).upper()==rr
 # Entire local alternate haplotype must be identical before and after normalization.
 old=lookup[r.id];start0=int(old['ROS_start0']);lo=min(start0,r.start)-50;hi=max(start0+len(old['ROS_ref']),r.start+len(rr))+50;reference=fa.fetch(r.chrom,lo,hi).upper();newhap=reference[:r.start-lo]+alt+reference[r.start-lo+len(rr):]
 candidates=[reference[:start0-lo]+a+reference[start0-lo+len(old['ROS_ref']):] for a in old['ROS_alt'].split(',')];assert newhap in candidates
 norm.append(dict(id=r.id,window_id=old['window_id'],chrom=r.chrom,start0=r.start,end0=r.start+len(rr),ref=rr,alt=alt,filter=old['filter'],source_record_max_alt_AF=old['max_alt_AF'],haplotype_equivalence_verified=True))
write('normalized_variant_alleles_ROS.tsv',norm)
# Include all annotated direct gene-child RNA features, not just one gene boundary.
features=[];genes={}
with (R.parent/'genomic.gff').open() as f:
 for l in f:
  if l.startswith('#'):continue
  v=l.rstrip().split('\t')
  if len(v)!=9:continue
  a={k:urllib.parse.unquote(x) for item in v[8].split(';') if '=' in item for k,x in [item.split('=',1)]}
  if v[2] in ['gene','pseudogene']:genes[a['ID']]=a.get('Name',a.get('gene',a['ID']))
  elif v[2] not in ['exon','CDS','region']:features.append((v,a))
tss=[];summary=[]
for w in rows:
 s,e=int(w['start']),int(w['end']);entries=[]
 for v,a in features:
  parents=a.get('Parent','').split(',')
  if v[0]!=w['chrom'] or not any(x in genes for x in parents) or v[6] not in ['+','-']:continue
  pos=int(v[3])-1 if v[6]=='+' else int(v[4])-1;gap=max(pos-e,s-(pos+1),0)
  if gap<=1000000:entries.append(dict(window_id=w['window_id'],chrom=v[0],feature=v[2],transcript_id=a.get('ID',''),gene=';'.join(genes[x] for x in parents if x in genes),strand=v[6],annotated_5prime0=pos,gap_bp=gap,window_overlap= s<=pos<e))
 entries.sort(key=lambda r:r['gap_bp']);tss+=entries;assert entries
 summary.append(dict(window_id=w['window_id'],nearest_annotated_transcript_5prime=entries[0],RNA_features_within_1Mb=len(entries),note='Annotated 5prime boundaries, not empirically mapped TSS; annotation incompleteness remains'))
write('annotated_transcript_5prime_context.tsv',tss)
write('variant_masks_PASS_ROS.tsv',[r for r in norm if r['filter']=='PASS'])
result=dict(status='completed',windows=2,projected_catalogue_records=len(projected),normalized_alleles=len(norm),all_REF_and_haplotypes_verified=True,transcript_5prime_summary=summary,bcftools_command=cmd,limitations=['Normalization is allele-level projection, not donor genotype validation','Maximum AF is source-record maximum, not individual split-allele AF','PASS variant REF span masks are not whole indel effect ranges','No guide or cut position selected','No chromatin contact data or functional regulatory validation'])
(O/'variant_regulatory_review.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
