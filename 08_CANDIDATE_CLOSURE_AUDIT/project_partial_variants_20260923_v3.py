from pathlib import Path
import json,pysam,subprocess,hashlib,datetime,collections
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';Q=C/'variant_catalogue_partial_query_20260923';O=C/'partial_variant_projection_20260923_v3'
s=json.loads((Q/'status.json').read_text());assert s['state']=='query_completed_projection_pending' and not s['errors']
O.mkdir(exist_ok=False);maps={r['id']:r for r in json.loads((C/'gapped_projection_20260923/UU_review.json').read_text())['results']}
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';U='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa';rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
projected={};pending=[];total=0
with pysam.FastaFile(str(R)) as ref,pysam.FastaFile(U) as uu:
 for item in s['completed']:
  d=json.loads((Q/(item['key']+'.json')).read_text());raw=Q/(item['key']+'.vcf');assert hashlib.sha256(raw.read_bytes()).hexdigest()==d['raw_vcf_sha256']
  w=maps[d['window']];assert w['mapq']>=30 and w['additional_reported_alignments']==0
  with pysam.VariantFile(str(raw)) as vf:rawa=[(r.pos,r.ref,alt,idx+1,None if r.info.get('AF') is None else r.info['AF'][idx]) for r in vf for idx,alt in enumerate(r.alts or ())]
  assert len(rawa)==len(d['alleles']);total+=len(rawa)
  for i,a in enumerate(d['alleles']):
   key=item['key']+'__'+str(i+1);assert rawa[i]==(a['pos1'],a['ref'],a['alt'],a['alt_index1'],a['AF']);p=a['pos1']-1;e=p+len(a['ref']);assert uu.fetch(w['uu_chrom'],p,e).upper()==a['ref'].upper()
   runs=[r for r in w['exact_runs'] if min(r['first_uu_pos0'],r['last_uu_pos0'])<=p and e<=max(r['first_uu_pos0'],r['last_uu_pos0'])+1]
   if not a['sequence_allele'] or len(runs)!=1:
    pending.append({'id':key,'allele':a,'reason':'nonsequence_ALT' if not a['sequence_allele'] else 'REF_not_inside_one_exact_run'});continue
   r=runs[0];rs=w['ros_start0']+r['ros_relative_start0'];re=w['ros_start0']+r['ros_relative_end0'];us=min(r['first_uu_pos0'],r['last_uu_pos0']);ue=max(r['first_uu_pos0'],r['last_uu_pos0'])+1
   u=uu.fetch(w['uu_chrom'],us,ue).upper();rr=ref.fetch(w['ros_chrom'],rs,re).upper();assert rr==(u if w['strand']=='+' else rc(u))
   pos=rs+p-us if w['strand']=='+' else re-(e-us);ar=a['ref'].upper();aa=a['alt'].upper()
   if w['strand']=='-':ar=rc(ar);aa=rc(aa)
   assert ref.fetch(w['ros_chrom'],pos,pos+len(ar)).upper()==ar
   uh=u[:p-us]+a['alt'].upper()+u[e-us:];rh=rr[:pos-rs]+aa+rr[pos-rs+len(ar):];assert rh==(uh if w['strand']=='+' else rc(uh))
   projected[key]={'window':w['id'],'chrom':w['ros_chrom'],'start0':pos,'ref':ar,'alt':aa,'AF':a['AF'],'filters':a['filters'],'exact_run_start0':rs,'exact_run_end0':re}
 head='##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c},length={ref.get_reference_length(c)}>\n' for c in sorted({x['chrom'] for x in projected.values()}))+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
 inp=O/'before_norm.vcf';out=O/'normalized.vcf';inp.write_text(head+''.join(f"{v['chrom']}\t{v['start0']+1}\t{k}\t{v['ref']}\t{v['alt']}\t.\t.\t.\n" for k,v in projected.items()))
 cmd=['/home/stardrako/miniforge3/envs/atac/bin/bcftools','norm','-f',str(R),'--check-ref','e','-m','-any','-Ov','-o',str(out),str(inp)];p=subprocess.run(cmd,capture_output=True,text=True);(O/'norm.log').write_text(p.stderr);assert p.returncode==0
 accepted=[];seen=set()
 with pysam.VariantFile(str(out)) as vf:
  for r in vf:
   assert r.id not in seen;seen.add(r.id);v=projected[r.id];p=v['start0'];assert ref.fetch(r.contig,r.start,r.stop).upper()==r.ref.upper()
   lo=max(0,min(p,r.start)-50);hi=min(ref.get_reference_length(r.contig),max(p+len(v['ref']),r.stop)+50)
   assert ref.fetch(r.contig,lo,p).upper()+v['alt']+ref.fetch(r.contig,p+len(v['ref']),hi).upper()==ref.fetch(r.contig,lo,r.start).upper()+r.alts[0].upper()+ref.fetch(r.contig,r.stop,hi).upper(), (v,str(r))
   if not(v['exact_run_start0']<=r.start and r.stop<=v['exact_run_end0']):pending.append({'id':r.id,'reason':'normalization_exits_exact_run'});continue
   accepted.append(dict(v,id=r.id,start0=r.start,ref=r.ref.upper(),alt=r.alts[0].upper()))
 assert seen==set(projected)
review={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_alleles':total,'accepted_partial_allele_records':len(accepted),'pending_records':len(pending),'pending_reasons':dict(collections.Counter(x['reason'] for x in pending)),'per_window':dict(collections.Counter(x['window'] for x in accepted)),'limitations':['Partial exact-run allele transfer only, never whole-window clearance','Source allele records not deduplicated across catalogues','Two multialignment windows not queried','Mismatches, indels and unmapped sequence remain unresolved for full-window variant assessment']}
(O/'review.json').write_text(json.dumps(review,indent=2));(O/'accepted_alleles.json').write_text(json.dumps(accepted,indent=2));(O/'pending.json').write_text(json.dumps(pending,indent=2));print(json.dumps(review,indent=2))


