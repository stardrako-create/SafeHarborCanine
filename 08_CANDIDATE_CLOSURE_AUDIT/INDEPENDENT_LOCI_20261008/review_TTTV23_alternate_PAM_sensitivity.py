from pathlib import Path
import json,gzip,pysam,hashlib
N=Path(__file__).resolve().parent;O=N/'hypothetical_TTTV23_panel_20261009';assert json.loads((O/'status.json').read_text())['state']=='completed';guides=json.loads((O/'guides_input.json').read_text());fa=pysam.FastaFile(str(N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];patterns=['TTTV','CTTV','TCTV','TTCV'];hits={k:{'ungapped':{},'bulge':{}} for k in guides}
def classify(pp):
 if pp is None:return 'unknown'
 if any(b not in 'ACGTN' for b in pp):return 'unknown'
 def matches(pattern):return all(b=='N' or b in ('ACG' if p=='V' else p) for b,p in zip(pp,pattern))
 if not any(matches(p) for p in patterns):return 'absent'
 return 'unknown' if 'N' in pp else ('canonical' if matches('TTTV') else 'single_T_to_C')
for mode in ['ungapped','bulge']:
 lookup=json.loads((O/(mode+'_query_map.json')).read_text());q={};k=None
 for line in (O/(mode+'_queries.fa')).read_text().splitlines():
  if line.startswith('>'):k=line[1:]
  else:q[k]=line
 with gzip.open(O/(mode+'_alignments.tsv.gz'),'rt') as f:
  for line in f:
   x=line.rstrip('\n').split('\t');gid=lookup[x[0]];g=guides[gid];c=x[2];s=int(x[3]);strand=x[1];length=len(q[x[0]])
   if mode=='ungapped' and c==g['chrom'] and s==g['start0'] and strand==g['strand']:continue
   if mode=='bulge' and c==g['chrom'] and s<g['start0']+23 and s+length>g['start0']:continue
   target=fa.fetch(c,s,s+length).upper();target=target if strand=='+' else rc(target);assert len(target)==length
   if mode=='ungapped':nm=sum(a!=b for a,b in zip(g['sequence'],target));assert nm<=3
   else:
    nm=min(sum(a!=b for a,b in zip(g['sequence'][:i]+g['sequence'][i+1:],target)) for i in range(23)) if length==22 else min(sum(a!=b for a,b in zip(g['sequence'],target[:i]+target[i+1:])) for i in range(24));assert nm<=2
   size=fa.get_reference_length(c);pp=fa.fetch(c,s-4,s).upper() if strand=='+' and s>=4 else (rc(fa.fetch(c,s+length,s+length+4).upper()) if strand=='-' and s+length+4<=size else None);cls=classify(pp)
   if cls=='absent':continue
   key=f'{c}:{s}:{length}:{strand}';old=hits[gid][mode].get(key)
   if old is None or nm<old['mismatches']:hits[gid][mode][key]=dict(chrom=c,start0=s,end0=s+length,strand=strand,PAM=pp,PAM_class=cls,mismatches=nm)
results=[]
for gid,g in guides.items():results.append(dict(guide_id=gid,window=g['window'],chrom=g['chrom'],start0=g['start0'],strand=g['strand'],sequence=g['sequence'],PAM=g['PAM'],ungapped_hits=list(hits[gid]['ungapped'].values()),single_bulge_hits=list(hits[gid]['bulge'].values())))
summary=[dict(window=w,exploratory_guides=sum(r['window']==w for r in results),no_hits_under_defined_reference_model=sum(r['window']==w and not r['ungapped_hits'] and not r['single_bulge_hits'] for r in results)) for w in sorted({r['window'] for r in results})]
out=N/'hypothetical_TTTV23_alternate_PAM_sensitivity_20261009.json';assert not out.exists();out.write_text(json.dumps(dict(state='completed',offtarget_PAM_patterns=patterns,source_primary='https://academic.oup.com/nar/article/48/7/3722/5763098',interpretation='Conservative diagnostic expansion of canonicalTTTV to oneT-to-C PAM substitutions documented forLbCas12aWT; not an exhaustive enzyme-specific PAM model or approved guides.',results=results,summary=summary,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [O/'results.json',O/'ungapped_alignments.tsv.gz',O/'bulge_alignments.tsv.gz']},limits=['No actual nuclease selected; orthologue/variant-dependent recognition and activity.','ROS reference only, ungapped<=3 substitutions or onebasebulge<=2 substitutions; other models/alternategenomes not covered.','Unknown PAMs retained; overlapping intendedintervalbulgealignments excluded explicitly.','Zero hits in this defined search is not zero biologicalofftargets or computationalauditcompletion.']),indent=2));print(json.dumps(summary))
