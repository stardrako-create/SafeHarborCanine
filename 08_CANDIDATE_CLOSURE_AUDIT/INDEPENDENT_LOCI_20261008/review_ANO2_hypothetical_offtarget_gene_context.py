from pathlib import Path
import json,gzip,pysam,hashlib,bisect,collections,urllib.parse
N=Path(__file__).resolve().parent;P=N.parent.parent;S=N/'hypothetical_NGG_panel_20261009/results.json';B=N/'hypothetical_NGG_single_bulge_20261009/results.json';K='NC_051831.1_39708074_39709074';panels=json.loads(S.read_text())['results'];bulges=json.loads(B.read_text())['results'];eligible={ (r['start0'],r['strand']):r for r in bulges if r['window']==K and sum(r['counts'][str(i)] for i in [0,1])==0};guides=[g for g in panels if g['window']==K and (g['start0'],g['strand']) in eligible];assert len(guides)==9
G=P/'01_referencia/ROS_Cfam_1.0/genomic.gff';bychrom=collections.defaultdict(list)
with G.open() as f:
 for line in f:
  if line.startswith('#'):continue
  v=line.rstrip('\n').split('\t')
  if len(v)!=9 or v[2] not in ['gene','pseudogene']:continue
  a={k:urllib.parse.unquote(x) for item in v[8].split(';') if '=' in item for k,x in [item.split('=',1)]};bychrom[v[0]].append((int(v[3])-1,int(v[4]),a.get('Name',a.get('ID','')),a.get('gene_biotype',''),v[6]))
indexes={}
for c,rr in bychrom.items():
 rr.sort();starts=[r[0] for r in rr];maxends=[];m=0
 for r in rr:m=max(m,r[1]);maxends.append(m)
 indexes[c]=(rr,starts,maxends)
def annotate(c,s,e):
 if c not in indexes:return []
 rr,starts,mends=indexes[c];i=bisect.bisect_left(starts,e)-1;hits=[]
 while i>=0 and mends[i]>s:
  r=rr[i]
  if r[1]>s and r[0]<e:hits.append(dict(gene=r[2],biotype=r[3],strand=r[4],gene_start0=r[0],gene_end0=r[1]))
  i-=1
 return hits
fa=pysam.FastaFile(str(P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];lookup=collections.defaultdict(list)
for g in guides:lookup[g['query']].append(g)
allhits={(g['start0'],g['strand']):[] for g in guides}
with gzip.open(N/'mismatch_specificity_20261009/alignments.tsv.gz','rt') as f:
 for line in f:
  x=line.rstrip('\n').split('\t')
  if x[0] not in lookup:continue
  c=x[2];s=int(x[3]);nm=0 if len(x)<8 or not x[7] else len(x[7].split(','))
  for g in lookup[x[0]]:
   if c==g['chrom'] and s==g['start0']:continue
   strand=x[1] if g['strand']=='+' else ('-' if x[1]=='+' else '+');size=fa.get_reference_length(c)
   pp=fa.fetch(c,s+20,s+23).upper() if strand=='+' and s+23<=size else (rc(fa.fetch(c,s-3,s).upper()) if strand=='-' and s>=3 else None)
   if pp is not None and (pp[1] not in 'GN' or pp[2] not in 'GN'):continue
   certainty='NGG' if pp is not None and set(pp)<=set('ACGT') else 'unknown_PAM';lo=s if strand=='+' else max(0,s-3);hi=min(size,s+23) if strand=='+' else s+20
   allhits[(g['start0'],g['strand'])].append(dict(chrom=c,start0=s,strand=strand,PAM=pp,PAM_status=certainty,mismatches=nm,gene_overlaps=annotate(c,lo,hi)))
results=[]
for g in guides:
 hh=allhits[(g['start0'],g['strand'])];br=eligible[(g['start0'],g['strand'])];bb=[dict(h,gene_overlaps=annotate(h['chrom'],h['start0'],h['end0'])) for h in br['one_bulge_NGG_or_unknown_hits']];assert len(hh)==sum(g['other_NGG_counts'].values())+sum(g['unknown_PAM_counts'].values());assert all(h['mismatches']==3 for h in hh)
 results.append(dict(guide=g,ungapped3mismatch_hits=hh,single_bulge2substitution_hits=bb,summary=dict(ungapped_hits=len(hh),ungapped_gene_overlapping_hits=sum(bool(h['gene_overlaps']) for h in hh),bulge_hits=len(bb),bulge_gene_overlapping_hits=sum(bool(h['gene_overlaps']) for h in bb))))
out=N/'ANO2_hypothetical_NGG_gene_context_20261009.json';assert not out.exists();out.write_text(json.dumps(dict(state='completed',guides_reviewed=9,selection_description='Exploratory NGG guides with no other NGG/unknown ungapped matches<=2 and no one-base-bulge matches<=1 substitutions; no new acceptance threshold imposed.',results=results,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S,B,G]},limits=['Gene overlap identifies annotation context, not editing propensity or gene harm.','All nine have other reference NGG3mismatch and onebulge2substitution matches.','No actual nuclease, alternate PAMs/genomes, empirical specificity, or biological clearance.']),indent=2));print(json.dumps([dict(start0=r['guide']['start0'],strand=r['guide']['strand'],**r['summary']) for r in results]))
