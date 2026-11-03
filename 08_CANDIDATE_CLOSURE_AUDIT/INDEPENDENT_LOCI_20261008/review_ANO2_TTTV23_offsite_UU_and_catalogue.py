from pathlib import Path
import json,pysam,subprocess,hashlib,collections,traceback,signal,datetime
N=Path(__file__).resolve().parent;S=N/'ANO2_TTTV23_one_SNP_sensitivity_20261009/results.json';O=N/'ANO2_TTTV23_offsite_UU_and_catalogue_20261009';O.mkdir(exist_ok=False);R=N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';U=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');ros=pysam.FastaFile(str(R));uu=pysam.FastaFile(str(U));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];sites={};inputs={};patterns={p+b for p in ['TTT','CTT','TCT','TTC'] for b in 'ACG'}
for g in json.loads(S.read_text())['results']:
 for h in g['one_SNP_sensitive_reference_sites']:
  key='site'+str(len(sites));c=h['chrom'];s=h['start0'];lo=max(0,s-500);hi=min(ros.get_reference_length(c),s+523);seq=ros.fetch(c,lo,hi).upper();sites[key]=dict(guide=g['guide'],site=h,ROS_context_start0=lo,ROS_context_end0=hi,context_sequence_sha256=hashlib.sha256(seq.encode()).hexdigest());inputs[key]=seq
assert len(sites)==23;(O/'contexts_ROS.fa').write_text(''.join('>'+k+'\n'+seq+'\n' for k,seq in inputs.items()));(O/'input_sites.json').write_text(json.dumps(sites,indent=2));status=dict(state='running',stage='mapping_to_UU',sites=len(sites),completed_queries=[],errors=[])
def save():(O/'status.json').write_text(json.dumps(status,indent=2))
def timeout(*a):raise TimeoutError('VCF range operation exceeded120seconds')
signal.signal(signal.SIGALRM,timeout);save();cmd=['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','4',str(U),str(O/'contexts_ROS.fa')];(O/'bwa.command.json').write_text(json.dumps(cmd))
with (O/'contexts_UU.sam').open('w') as out,(O/'bwa.err').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True,timeout=1200)
records=collections.defaultdict(list)
with pysam.AlignmentFile(str(O/'contexts_UU.sam'),'r') as sam:
 for r in sam.fetch(until_eof=True):records[r.query_name].append(r)
results=[]
for key,j in sites.items():
 rr=records[key];primary=[r for r in rr if not r.is_secondary and not r.is_supplementary and not r.is_unmapped];result=dict(site_id=key,**j,projection_status='not_resolved')
 if len(primary)!=1:result['projection_reason']='No unique primary';results.append(result);continue
 r=primary[0];alternatives=[a for a in rr if not a.is_unmapped and (a.is_secondary or a.is_supplementary)];result['UU_alignment']=dict(chrom=r.reference_name,start0=r.reference_start,end0=r.reference_end,cigar=r.cigarstring,mapq=r.mapping_quality,NM=r.get_tag('NM') if r.has_tag('NM') else None,reverse=r.is_reverse,alternative_count=len(alternatives))
 if r.mapping_quality<30 or alternatives or r.has_tag('XA') or r.has_tag('SA'):result['projection_reason']='Mapping rule not satisfied';results.append(result);continue
 mapping={}
 for qi,ri in r.get_aligned_pairs(matches_only=False):
  if qi is None or ri is None:continue
  orig=qi if not r.is_reverse else len(inputs[key])-1-qi
  assert orig not in mapping;mapping[orig]=ri
 h=j['site'];core_lo=h['start0']-4 if h['strand']=='+' else h['start0'];core_hi=h['end0'] if h['strand']=='+' else h['end0']+4;offset=core_lo-j['ROS_context_start0'];targetpos=[mapping.get(i) for i in range(offset,offset+27)]
 if core_lo<0 or core_hi>ros.get_reference_length(h['chrom']) or any(p is None for p in targetpos) or len(set(targetpos))!=27:result['projection_reason']='Target27bp not fully mapped';results.append(result);continue
 step=-1 if r.is_reverse else 1
 if any(b-a!=step for a,b in zip(targetpos,targetpos[1:])):result['projection_reason']='Gapped/noncontiguous target27bp';results.append(result);continue
 lo=min(targetpos);hi=max(targetpos)+1;strand=h['strand'] if not r.is_reverse else ('-' if h['strand']=='+' else '+');dna=uu.fetch(r.reference_name,lo,hi).upper();oriented=dna if strand=='+' else rc(dna);assert len(oriented)==27;nm=sum(a!=b for a,b in zip(j['guide']['protospacer_5to3'],oriented[4:]));modelhit=oriented[:4] in patterns and nm<=3
 result.update(projection_status='unique_contiguous_target27bp',UU_core=dict(chrom=r.reference_name,start0=lo,end0=hi,target_strand=strand,sequence_plus=dna,sequence_oriented=oriented,UU_PAM=oriented[:4],UU_guide_mismatches=nm,UU_reference_hit_in_defined_model=modelhit),catalogue=[]);results.append(result)
(O/'projection_review.json').write_text(json.dumps(dict(results=results,source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),sam_sha256=hashlib.sha256((O/'contexts_UU.sam').read_bytes()).hexdigest(),limits=['Unique BWAcontext mapping plus contiguous27baseprojection; nonidenticalbases not silentlytreatedidentical.','UUreference comparison is another assembly sequence, not empirical cleavage or independentdogvalidation.']),indent=2));status['stage']='catalogue_range_queries';save();base='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/'
for src in ['AutoAndXPAR.SNPs.vqsr99.vcf.gz','AutoAndXPAR.nonSNPs.filter.vcf.gz']:
 try:signal.alarm(120);vf=pysam.VariantFile(base+src,drop_samples=True);signal.alarm(0);(O/(src+'.header.txt')).write_text(str(vf.header))
 except Exception as e:signal.alarm(0);status['errors'].append(dict(source=src,error=str(e)));save();continue
 for result in results:
  if result['projection_status']!='unique_contiguous_target27bp':continue
  key=result['site_id'];core=result['UU_core'];aliases=[c for c in vf.header.contigs if vf.header.contigs[c].length==uu.get_reference_length(core['chrom'])]
  if len(aliases)!=1:result['catalogue'].append(dict(source=src,status='unavailable_alias',aliases=aliases));continue
  alias=aliases[0]
  if alias.lower().replace('chr','') in ['x','y','m','mt']:result['catalogue'].append(dict(source=src,status='not_covered_by_AutoAndXPAR_for_this_nonautosomal_site',alias=alias));continue
  try:
   signal.alarm(120);rr=list(vf.fetch(alias,core['start0'],core['end0']));signal.alarm(0);raw=O/(key+'__'+src+'.vcf');raw.write_text(str(vf.header)+''.join(str(r) for r in rr));alleles=[]
   for r in rr:
    observed=uu.fetch(core['chrom'],r.start,r.start+len(r.ref)).upper();assert observed==r.ref.upper(),'CatalogueREF differsfromUUreference';af=r.info.get('AF');ac=r.info.get('AC');an=r.info.get('AN')
    if af is not None and not isinstance(af,(list,tuple)):af=(af,)
    if ac is not None and not isinstance(ac,(list,tuple)):ac=(ac,)
    if af is not None:assert len(af)==len(r.alts or ())
    if ac is not None:assert len(ac)==len(r.alts or ())
    for ai,alt in enumerate(r.alts or ()):
     a=dict(pos1=r.pos,ref=r.ref,alt=alt,filters=list(r.filter),AF=None if af is None else af[ai],AC=None if ac is None else ac[ai],AN=an,REF_matches_UU=True,status='nonSNP_or_boundary_effect_not_assessed')
     if len(r.ref)==1 and len(alt)==1 and set(alt.upper())<=set('ACGT') and core['start0']<=r.start<core['end0']:
      off=r.start-core['start0'];dna=core['sequence_plus'];assert dna[off]==r.ref.upper();mutant=dna[:off]+alt.upper()+dna[off+1:];oriented=mutant if core['target_strand']=='+' else rc(mutant);nm=sum(x!=y for x,y in zip(result['guide']['protospacer_5to3'],oriented[4:]));a.update(status='single_SNP_model_assessed',alternate_PAM=oriented[:4],alternate_guide_mismatches=nm,alternate_hit_in_defined_model=oriented[:4] in patterns and nm<=3)
     alleles.append(a)
   query=dict(source=base+src,status='range_query_completed',alias=alias,record_count=len(rr),alleles=alleles,raw_vcf_sha256=hashlib.sha256(raw.read_bytes()).hexdigest());result['catalogue'].append(query);status['completed_queries'].append(dict(site_id=key,source=src,records=len(rr)));save()
  except Exception as e:signal.alarm(0);result['catalogue'].append(dict(source=src,status='failed',error=str(e)));status['errors'].append(dict(site_id=key,source=src,error=str(e)));save()
 vf.close()
summary=[]
for g in json.loads(S.read_text())['results']:
 name=g['guide']['hypothesis_id'];subset=[r for r in results if r['guide']['hypothesis_id']==name];summary.append(dict(guide=name,sensitive_sites=len(subset),unique_contiguous_UU_targets=sum(r['projection_status']=='unique_contiguous_target27bp' for r in subset),UU_reference_hits=sum(r.get('UU_core',{}).get('UU_reference_hit_in_defined_model',False) for r in subset),PASS_single_SNP_model_hits=sum(a.get('alternate_hit_in_defined_model',False) and a['filters']==['PASS'] for r in subset for c in r.get('catalogue',[]) for a in c.get('alleles',[])),nonPASS_single_SNP_model_hits=sum(a.get('alternate_hit_in_defined_model',False) and a['filters']!=['PASS'] for r in subset for c in r.get('catalogue',[]) for a in c.get('alleles',[])),unassessed_nonSNP_or_boundary_alleles=sum(a['status']!='single_SNP_model_assessed' for r in subset for c in r.get('catalogue',[]) for a in c.get('alleles',[]))))
result=dict(results=results,summary=summary,source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),limits=['SNPs applied singly toUUreference background; no observedhaplotypephase orcombinedvariantclearance.','Unmapped/nonautosomal/nonSNP/boundarycases remainexplicitlyunassessed, neverPASSfrommissingdata.','This is reference/catalogueoffsite sensitivity, not experimental cleavage or finalguideapproval.','No complete alternate-genome/populationSVorbulgeplusvariantsearch.']);(O/'results.json').write_text(json.dumps(result,indent=2));status.update(state='completed' if not status['errors'] else 'completed_with_errors',stage='review_complete',summary=summary,utc_end=datetime.datetime.now(datetime.timezone.utc).isoformat());save();print(json.dumps(status))
