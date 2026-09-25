from pathlib import Path
import csv,json,hashlib,datetime,subprocess,collections
import pysam
A=Path(__file__).resolve().parent; C=A/'CONTROL_LOCUS_REVIEW_20260922'; Q=C/'variant_catalogue_query_20260923'; O=C/'variant_projection_20260923';O.mkdir(exist_ok=False)
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'
U='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa'
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
rows={r['id']:r for r in csv.DictReader((C/'transcript_and_boundary_review.tsv').open(),delimiter='\t')}
status=json.loads((Q/'status.json').read_text());assert status['state']=='query_completed_projection_pending' and len(status['completed'])==16 and not status['errors']
orig={};pending=[];sources=[];totalrecords=0
with pysam.FastaFile(str(R)) as ref,pysam.FastaFile(U) as uu:
 for item in status['completed']:
  f=Q/(item['key']+'.json');d=json.loads(f.read_text());raw=Q/(item['key']+'.vcf')
  assert hashlib.sha256(raw.read_bytes()).hexdigest()==d['raw_vcf_sha256']
  with pysam.VariantFile(str(raw)) as vf:
   records=list(vf)
  assert len(records)==d['record_count']==item['records'];totalrecords+=len(records)
  rawalleles=[(r,idx,alt) for r in records for idx,alt in enumerate(r.alts or ())]
  assert len(rawalleles)==len(d['alleles'])==item['alleles']
  h=d['mapping'];w=rows[d['window']];rs=int(w['start0']);re=int(w['end0']);chrom=w['chrom']
  useq=uu.fetch(h['chrom'],h['start0'],h['end0']).upper();rseq=ref.fetch(chrom,rs,re).upper()
  assert rseq==(useq if h['strand']=='+' else rc(useq)) and len(rseq)==1000
  for n,(r,idx,alt) in enumerate(rawalleles):
   a=d['alleles'][n];af=r.info.get('AF');af=None if af is None else af[idx]
   assert (r.pos,r.ref,alt,idx+1,af)==(a['pos1'],a['ref'],a['alt'],a['alt_index1'],a['AF'])
   assert uu.fetch(h['chrom'],r.start,r.start+len(r.ref)).upper()==r.ref.upper()
   key=item['key']+'__'+str(n+1)
   if not a['within_exact_window'] or not a['sequence_allele'] or set(r.ref.upper())-set('ACGT'):
    pending.append({'id':key,'allele':a,'reason':'outside exact window or noncanonical allele'});continue
   p=rs+r.start-h['start0'] if h['strand']=='+' else re-(r.start-h['start0']+len(r.ref))
   rr=r.ref.upper() if h['strand']=='+' else rc(r.ref.upper());aa=alt.upper() if h['strand']=='+' else rc(alt.upper())
   assert ref.fetch(chrom,p,p+len(rr)).upper()==rr
   ua=useq[:r.start-h['start0']]+alt.upper()+useq[r.start-h['start0']+len(r.ref):]
   ra=rseq[:p-rs]+aa+rseq[p-rs+len(rr):]
   assert ra==(ua if h['strand']=='+' else rc(ua))
   orig[key]={'window':d['window'],'chrom':chrom,'start0':p,'ref':rr,'alt':aa,'AF':af,'filters':list(r.filter),'source_key':item['key']}
  sources.append({'file':f.name,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
 header='##fileformat=VCFv4.2\n'+''.join(f'##contig=<ID={c},length={ref.get_reference_length(c)}>\n' for c in sorted({v['chrom'] for v in orig.values()}))+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
 inp=O/'projected_before_norm.vcf';out=O/'projected_normalized.vcf'
 inp.write_text(header+''.join(f"{v['chrom']}\t{v['start0']+1}\t{k}\t{v['ref']}\t{v['alt']}\t.\t.\t.\n" for k,v in orig.items()))
 cmd=['/home/stardrako/miniforge3/envs/atac/bin/bcftools','norm','-f',str(R),'--check-ref','e','-m','-any','-Ov','-o',str(out),str(inp)]
 proc=subprocess.run(cmd,capture_output=True,text=True);(O/'bcftools.log').write_text(proc.stderr);assert proc.returncode==0,proc.stderr
 normalized=[];seen=set()
 with pysam.VariantFile(str(out)) as vf:
  for r in vf:
   assert r.id not in seen;seen.add(r.id);v=orig[r.id];p=v['start0'];assert r.contig==v['chrom'] and len(r.alts)==1
   assert ref.fetch(r.contig,r.start,r.stop).upper()==r.ref.upper()
   lo=max(0,min(p,r.start)-50);hi=min(ref.get_reference_length(r.contig),max(p+len(v['ref']),r.stop)+50)
   before=ref.fetch(r.contig,lo,p).upper()+v['alt']+ref.fetch(r.contig,p+len(v['ref']),hi).upper()
   after=ref.fetch(r.contig,lo,r.start).upper()+r.alts[0]+ref.fetch(r.contig,r.stop,hi).upper();assert before==after
   w=rows[v['window']]
   normalized.append(dict(v,id=r.id,start0=r.start,ref=r.ref,alt=r.alts[0],normalized_within_window=r.start>=int(w['start0']) and r.stop<=int(w['end0'])))
 assert seen==set(orig)
summary={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'state':'completed','source_records':totalrecords,'source_alleles':sum(x['alleles'] for x in status['completed']),'projected_normalized_allele_records':len(normalized),'pending_alleles':pending,'source_receipts':sources,'per_window':{},'checks':['Raw VCF hashes and counts','AF checked per ALT against raw VCF','UU and ROS REF exact','Whole-window alternate haplotype equivalent across orientation','Normalized REF and flanking alternate haplotype equivalent','One normalized record per source allele ID'],'limitations':['Catalogue alleles, not donor genotypes','Counts retain source allele records; not deduplicated across catalogues','12 other windows remain unresolved','No guide or biological validation']}
for w in sorted({a['window'] for a in normalized}):
 arr=[a for a in normalized if a['window']==w];summary['per_window'][w]={'allele_records':len(arr),'PASS':sum(a['filters']==['PASS'] for a in arr),'PASS_AF_ge_001':sum(a['filters']==['PASS'] and a['AF'] is not None and a['AF']>=.01 for a in arr),'normalized_outside_window':sum(not a['normalized_within_window'] for a in arr)}
(O/'alleles.json').write_text(json.dumps(normalized,indent=2));(O/'review.json').write_text(json.dumps(summary,indent=2));print(json.dumps({k:v for k,v in summary.items() if k not in ['source_receipts','pending_alleles']},indent=2))
