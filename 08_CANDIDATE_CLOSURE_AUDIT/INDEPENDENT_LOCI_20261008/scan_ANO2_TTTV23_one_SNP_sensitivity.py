from pathlib import Path
import json,pysam,subprocess,gzip,hashlib,collections,traceback
N=Path(__file__).resolve().parent;S=N/'ANO2_TTTV23_REFERENCE_MODEL_CANDIDATES_20261009.json';O=N/'ANO2_TTTV23_one_SNP_sensitivity_20261009';O.mkdir(exist_ok=False);guides=json.loads(S.read_text())['results'];assert len(guides)==2;fa=pysam.FastaFile(str(N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];queries={};lookup={};patterns=[prefix+b for prefix in ['TTT','CTT','TCT','TTC'] for b in 'ACG']
for i,g in enumerate(guides):
 seq=g['protospacer_5to3']
 for pos in range(23):
  for b in 'ACGT':
   if seq[pos]==b:continue
   key='q'+str(len(queries));queries[key]=seq[:pos]+b+seq[pos+1:];lookup[key]=i
(O/'queries.fa').write_text(''.join('>'+k+'\n'+v+'\n' for k,v in queries.items()));(O/'query_map.json').write_text(json.dumps(lookup));cmd=['/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin/bowtie','-f','-v','3','-a','--threads','4','/mnt/d/Jin2024_work/bowtie1_ROS_20261009/ROS',str(O/'queries.fa')];(O/'command.json').write_text(json.dumps(cmd));status=dict(state='running',queries=len(queries),reference_model='Enumerate all69single-substitutiontemplatesper23ntguide withBowtie<=3; complete ACGTreference search up to4 originalguide mismatches; canonical/oneT>C PAM sensitivity',raw_alignments=0);(O/'status.json').write_text(json.dumps(status));hits={i:{} for i in range(2)}
try:
 with gzip.open(O/'raw_alignments.tsv.gz','wt',compresslevel=1) as raw,(O/'scan.err').open('w') as err:
  p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True)
  for line in p.stdout:
   raw.write(line);x=line.rstrip('\n').split('\t');idx=lookup[x[0]];g=guides[idx];c=x[2];s=int(x[3]);strand=x[1];target=fa.fetch(c,s,s+23).upper();target=target if strand=='+' else rc(target);assert len(target)==23 and sum(a!=b for a,b in zip(queries[x[0]],target))<=3;nm=sum(a!=b for a,b in zip(g['protospacer_5to3'],target));assert nm<=4;status['raw_alignments']+=1
   if c==g['chrom'] and s==g['protospacer_start0'] and strand==g['strand']:continue
   size=fa.get_reference_length(c);pp=fa.fetch(c,s-4,s).upper() if strand=='+' and s>=4 else (rc(fa.fetch(c,s+23,s+27).upper()) if strand=='-' and s+27<=size else None)
   if pp is None or set(pp)-set('ACGT'):dist=None
   else:dist=min(sum(a!=b for a,b in zip(pp,pat)) for pat in patterns)
   relevant=(nm==4 and dist==0) or (nm<=3 and dist==1) or (nm<=3 and dist is None)
   if not relevant:continue
   key=f'{c}:{s}:{strand}';hits[idx][key]=dict(chrom=c,start0=s,end0=s+23,strand=strand,oriented_protospacer=target,PAM=pp,original_guide_mismatches=nm,PAM_distance_to_model=dist,sensitivity_case='one_protospacer_SNP_could_reduce4to3' if nm==4 else ('one_PAM_SNP_could_create_recognized_family' if dist==1 else 'unknown_PAM'),mismatch_positions=[i for i,(a,b) in enumerate(zip(g['protospacer_5to3'],target)) if a!=b])
  assert p.wait()==0
 result=[dict(guide=g,one_SNP_sensitive_reference_sites=list(hits[i].values()),sites=len(hits[i])) for i,g in enumerate(guides)];(O/'results.json').write_text(json.dumps(dict(results=result,source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),limits=['Sensitivity to a hypothetical singleSNP, not evidence that such alleles exist in dogs or cause cleavage.','No populationVCF/genotype verification in this stage.','No full alternate-genome, structuralvariation or bulge+SNP enumeration.','No actual nuclease/guideapproval; PAMmodel limited toexplicitfamily.','AssemblyN/gaps aremissinginformation.']),indent=2));status.update(state='completed',summary=[dict(guide=g['hypothesis_id'],one_SNP_sensitive_sites=len(hits[i])) for i,g in enumerate(guides)])
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:(O/'status.json').write_text(json.dumps(status,indent=2));print(json.dumps(status))
