from pathlib import Path
import json,csv,pysam,hashlib,datetime
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=A/'W01_FOCUSED_REVIEW_20260923';O.mkdir(exist_ok=False)
ids={'bg1k_2848','bg1k_0442'};integrated=json.loads((C/'integrated_review_20260923/review.json').read_text())
for f,h in integrated['sources_sha256'].items():assert hashlib.sha256((C/f).read_bytes()).hexdigest()==h
rows={r['id']:r for r in csv.DictReader((C/'transcript_and_boundary_review.tsv').open(),delimiter='\t') if r['id'] in ids};alleles=json.loads((C/'variant_projection_20260923/alleles.json').read_text());segments=list(csv.DictReader((C/'exact_sequences/all_short_segments_exact_uniqueness.tsv').open(),delimiter='\t'));seqreview={}
for key in ids:
 variants=[a for a in alleles if a['window']==key and a['filters']==['PASS']]
 seg=[r for r in segments if r['window_id']==key and int(r['length'])==20];assert len(seg)==981
 unique=[r for r in seg if r['unique_exact']=='True'];clear=[r for r in unique if not any(int(r['start0'])<a['start0']+len(a['ref']) and a['start0']<int(r['end0']) for a in variants)]
 seqreview[key]={'exact_unique_20mers':len(unique),'exact_unique_PASS_REF_clear_20mers':len(clear),'PASS_allele_records':len(variants),'scope':'Exact-locus count and catalogued REF span screen only; no PAM, mismatches, bulges, individual genotype or validated guides'}
svfile=Path('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz');maps=json.loads((C/'UU_projection/projection_review.json').read_text())['results'];queries=[]
for w in maps:
 if w['id'] in ids:
  h=next(h for h in w['hits'] if h['primary'] and h['exact_sequence_verified']);queries.append(dict(h,id=w['id']))
queries.append({'id':'w01','chrom':'NC_049228.1','start0':17403010,'end0':17404010})
svresults=[]
with pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa') as ref,pysam.VariantFile(str(svfile)) as vf:
 for q in queries:
  aliases=[c for c in vf.header.contigs if vf.header.contigs[c].length==ref.get_reference_length(q['chrom'])];assert len(aliases)==1
  found=[]
  for r in vf.fetch(aliases[0],q['start0'],q['end0']):
   if r.start<q['end0'] and r.stop>q['start0']:found.append({'pos1':r.pos,'end0':r.stop,'id':r.id,'ref':r.ref,'alts':r.alts,'filter':list(r.filter),'SVTYPE':r.info.get('SVTYPE')})
  svresults.append({'id':q['id'],'uu_contig':q['chrom'],'vcf_contig':aliases[0],'start0':q['start0'],'end0':q['end0'],'overlap_records':found,'scope':'Indexed interval-overlap query; not BND mate or breakpoint confidence-interval analysis; catalogue absence does not exclude individual SV'})
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'controls':rows,'sequence_screen':seqreview,'SV_overlap_query':svresults,'SV_source':str(svfile),'SV_file_size':svfile.stat().st_size,'integrated_source_hashes_verified':len(integrated['sources_sha256'])}
(O/'focused_checks.json').write_text(json.dumps(res,indent=2));print(json.dumps({'sequence':seqreview,'SV':svresults},indent=2))
