from pathlib import Path
import json,gzip,pysam,subprocess,hashlib,collections,traceback,fcntl
N=Path(__file__).resolve().parent;O=N/'hypothetical_TTTV23_panel_20261009';O.mkdir(exist_ok=False);lock=open('/tmp/safeharbor_TTTV23_20261009.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
R=N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';fa=pysam.FastaFile(str(R));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];A=N.parent/'ALL_CANDIDATES_20261005';masks={r['window']:set(r['mask_offsets']) for r in json.loads((N/'two_loci_REF_mask_specificity.json').read_text())['results']};ano='NC_051831.1_39708074_39709074';masks[ano]=set(next(r for r in json.loads((A/'masked_specificity_20261008.json').read_text())['results'] if r['id']==ano)['masked_offsets'])
def pam(c,s,l,strand):
 size=fa.get_reference_length(c)
 return fa.fetch(c,s-4,s).upper() if strand=='+' and s>=4 else (rc(fa.fetch(c,s+l,s+l+4).upper()) if strand=='-' and s+l+4<=size else None)
def classify(pp):
 if pp is None or any(b not in 'ACGTN' for b in pp):return 'unknown'
 if any(b not in 'TN' for b in pp[:3]) or pp[3] not in 'ACGN':return 'absent'
 return 'TTTV' if 'N' not in pp else 'unknown'
guides={}
for w,mask in masks.items():
 c,start,end=w.rsplit('_',2);start=int(start);end=int(end)
 for s in range(start,end-22):
  for strand in '+-':
   lo=s-4 if strand=='+' else s;hi=s+23 if strand=='+' else s+27
   if lo<start or hi>end or any(i in mask for i in range(lo-start,hi-start)):continue
   pp=pam(c,s,23,strand)
   if classify(pp)!='TTTV':continue
   seq=fa.fetch(c,s,s+23).upper();seq=seq if strand=='+' else rc(seq)
   if set(seq)-set('ACGT'):continue
   guides['g'+str(len(guides))]=dict(window=w,chrom=c,start0=s,strand=strand,sequence=seq,PAM=pp,exact_self_hits=0,ungapped_hits={},single_bulge_hits={})
(O/'guides_input.json').write_text(json.dumps(guides,indent=2));status=dict(state='running',hypothesis='Exploratory23nt target with5prime canonicalTTTV; not actual nuclease selection',guides=len(guides),raw_alignments=0,source_urls=['https://www.neb.com/en-us/products/m0653-engen-lba-cas12a-cpf1','https://www.neb.com/faqs/how-do-i-design-a-guide-rna-for-use-with-engen-lba-cas12a'])
def save():(O/'status.json').write_text(json.dumps(status,indent=2))
def scan(mode):
 queries={};lookup={}
 for gid,g in guides.items():
  seq=g['sequence'];templates={seq} if mode=='ungapped' else set(seq[:i]+seq[i+1:] for i in range(23))|set(seq[:i]+b+seq[i:] for i in range(24) for b in 'ACGT')
  for t in sorted(templates):key='q'+str(len(queries));queries[key]=t;lookup[key]=gid
 qf=O/(mode+'_queries.fa');qf.write_text(''.join('>'+k+'\n'+v+'\n' for k,v in queries.items()));(O/(mode+'_query_map.json')).write_text(json.dumps(lookup));cmd=['/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin/bowtie','-f','-v','3' if mode=='ungapped' else '2','-a','--threads','4','/mnt/d/Jin2024_work/bowtie1_ROS_20261009/ROS',str(qf)];(O/(mode+'.command.json')).write_text(json.dumps(cmd));status.update(mode=mode,queries=len(queries));save()
 with gzip.open(O/(mode+'_alignments.tsv.gz'),'wt',compresslevel=1) as raw,(O/(mode+'.err')).open('w') as err:
  p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=err,text=True)
  for line in p.stdout:
   raw.write(line);x=line.rstrip('\n').split('\t');g=guides[lookup[x[0]]];c=x[2];s=int(x[3]);strand=x[1];length=len(queries[x[0]]);target=fa.fetch(c,s,s+length).upper();target=target if strand=='+' else rc(target);nm=sum(a!=b for a,b in zip(queries[x[0]],target));assert len(target)==length and nm<=(3 if mode=='ungapped' else 2);status['raw_alignments']+=1
   if mode=='ungapped' and c==g['chrom'] and s==g['start0'] and strand==g['strand']:assert nm==0;g['exact_self_hits']+=1;continue
   if mode=='bulge':
    nm=min(sum(a!=b for a,b in zip(g['sequence'][:i]+g['sequence'][i+1:],target)) for i in range(23)) if length==22 else min(sum(a!=b for a,b in zip(g['sequence'],target[:i]+target[i+1:])) for i in range(24));assert nm<=2
    if c==g['chrom'] and s<g['start0']+23 and s+length>g['start0']:continue
   pp=pam(c,s,length,strand);cls=classify(pp)
   if cls=='absent':continue
   field='ungapped_hits' if mode=='ungapped' else 'single_bulge_hits';key=f'{c}:{s}:{length}:{strand}';old=g[field].get(key)
   if old is None or nm<old['mismatches']:g[field][key]=dict(chrom=c,start0=s,end0=s+length,strand=strand,PAM=pp,PAM_status=cls,mismatches=nm,bulge=('RNA' if length==22 else 'DNA') if mode=='bulge' else None)
   if status['raw_alignments']%50000==0:save()
  assert p.wait()==0
try:
 scan('ungapped');assert all(g['exact_self_hits']==1 for g in guides.values());scan('bulge');results=[]
 for gid,g in guides.items():
  g=dict(g,guide_id=gid,ungapped_hits=list(g['ungapped_hits'].values()),single_bulge_hits=list(g['single_bulge_hits'].values()));results.append(g)
 summaries=[]
 for w in sorted(masks):
  rr=[r for r in results if r['window']==w];summaries.append(dict(window=w,guides=len(rr),no_ungapped_hits_through3_and_no_single_bulge_hits_through2=sum(not r['ungapped_hits'] and not r['single_bulge_hits'] for r in rr),no_ungapped_hits_through3=sum(not r['ungapped_hits'] for r in rr)))
 (O/'results.json').write_text(json.dumps(dict(results=results,summary=summaries,limits=['Hypothetical canonicalTTTV23nt panel only, not actual nuclease or approved guides.','No alternate PAMrecognition, double/multibasebulges or populationalternate-genomeofftargets.','Bulge hits overlapping intended23ntinterval excluded; not local shiftedcleavageclearance.','Not biologicalvalidation or new locus acceptance threshold.']),indent=2));status.update(state='completed',summary=summaries)
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc())
finally:save();print(json.dumps(status))
