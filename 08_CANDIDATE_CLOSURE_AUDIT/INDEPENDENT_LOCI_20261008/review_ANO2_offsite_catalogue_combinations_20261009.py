from pathlib import Path
import json,pysam,hashlib,itertools,datetime
N=Path(__file__).resolve().parent;source=N/'review_ANO2_offsite_three_indels_20261009.py';text=source.read_text(encoding='utf-8-sig');functions=text[text.index('rc=lambda'):text.index('out=[]')];exec(functions)
S=N/'ANO2_TTTV23_offsite_UU_and_catalogue_20261009/results.json';O=N/'ANO2_offsite_catalogue_combinations_20261009';O.mkdir(exist_ok=False);f=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');results=[]
for r in json.loads(S.read_text())['results']:
 core=r['UU_core'];alleles=[a for q in r['catalogue'] for a in q.get('alleles',[])];assert len(alleles)<=3
 lo=max(0,core['start0']-100);hi=core['end0']+100;ref=f.fetch(core['chrom'],lo,hi).upper();synthetic=[]
 for flags in itertools.product((0,1),repeat=len(alleles)):
  chosen=sorted([a for a,flag in zip(alleles,flags) if flag],key=lambda a:a['pos1']);assert all(a['pos1']-1+len(a['ref'])<=b['pos1']-1 for a,b in zip(chosen,chosen[1:]))
  seq=ref
  for a in reversed(chosen):
   offset=a['pos1']-1-lo;assert seq[offset:offset+len(a['ref'])]==a['ref'];seq=seq[:offset]+a['alt']+seq[offset+len(a['ref']):]
  hits=scan(seq,r['guide']['protospacer_5to3']);synthetic.append(dict(selected_alleles=chosen,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),model_hits=hits))
 results.append(dict(site=r['site_id'],guide=r['guide']['hypothesis_id'],chrom=core['chrom'],start0=lo,end0=hi,alleles=len(alleles),combinations=synthetic))
receipt=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),scan_function_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),results=results,synthetic_contexts=sum(len(r['combinations']) for r in results),contexts_with_model_hits=sum(bool(c['model_hits']) for r in results for c in r['combinations']),limits=['All combinations of only the alleles retrieved in the23 selected27bp sites; variants elsewhere and flanking alleles not queried.','Synthetic combinations, not observed/phased canine haplotypes. PASS and nonPASS provenance retained, no nonPASS reinterpretation as verified allele.','Defined exploratory PAM model, not exhaustive specificity or empirical validation.'])
(O/'receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({k:receipt[k] for k in ['state','synthetic_contexts','contexts_with_model_hits']}))
