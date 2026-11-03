from pathlib import Path
import csv,json,gzip,urllib.parse,hashlib,datetime,collections
N=Path(__file__).resolve().parent;P=N.parent.parent;riskfile=P/'05_SHIP/canine_risk_genes.tsv';risk={r['gene_symbol']:r['categories'] for r in csv.DictReader(riskfile.open(),delimiter='\t')};R=P/'01_referencia/ROS_Cfam_1.0/genomic.gff';U=N/'ANO2_UU_annotation_crosscheck_20261009/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gff.gz'
def genes(path):
 opener=gzip.open if path.suffix=='.gz' else open
 with opener(path,'rt') as f:
  for line in f:
   if line.startswith('#'):continue
   x=line.rstrip().split('\t')
   if len(x)!=9 or x[2]!='gene':continue
   a={k:urllib.parse.unquote(v) for item in x[8].split(';') if '=' in item for k,v in [item.split('=',1)]};names={a.get('gene',''),a.get('Name','')}|set(a.get('gene_synonym','').split(','));ids={item.split(':',1)[1] for item in a.get('Dbxref','').split(',') if item.startswith('GeneID:')};yield x,a,names-{''},ids
byid=collections.defaultdict(set);mappedROS=set();direct=[]
for x,a,names,ids in genes(R):
 hit=names&set(risk)
 if hit:mappedROS|=hit
 for gid in ids:byid[gid]|=hit
mappedUU=set();nearest=[];lo=6741835;hi=6742835
for x,a,names,ids in genes(U):
 direct=names&set(risk);via=set().union(*(byid.get(gid,set()) for gid in ids));matched=direct|via;mappedUU|=matched
 if x[0]!='NC_049248.1' or not matched:continue
 start=int(x[3])-1;end=int(x[4]);d=lo-end+1 if end<=lo else (start-hi+1 if start>=hi else 0)
 nearest.append(dict(UU_gene=a.get('gene'),attributes=a,matched_risk_symbols=sorted(matched),direct_or_alias_symbols=sorted(direct),ROS_GeneID_linked_symbols=sorted(via),geneids=sorted(ids),start0=start,end0=end,distance_bp=d))
assert nearest
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),risk_symbols=len(risk),mapped_via_ROS_symbols_or_aliases=len(mappedROS),mapped_UU_with_alias_and_GeneID=len(mappedUU),unresolved_symbols=sorted(set(risk)-mappedUU),nearest_UU_risk_genes=sorted(nearest,key=lambda r:r['distance_bp'])[:20],risk_hits_within300kb=[r for r in nearest if r['distance_bp']<300000],sources_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [riskfile,R,U]},script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limits=['Frozen risk list, exact symbols/annotated synonyms plus stable GeneID links from ROS; no inferred orthologue identity or risk-list expansion.','Unresolved symbols are missing mapping evidence, not confirmed absence of risk.','Another assembly/annotation view, not independent biological confirmation.'])
p=N/'ANO2_UU_risk_alias_GeneID_crosscheck_20261009.json';assert not p.exists();p.write_text(json.dumps(result,indent=2));print(json.dumps(dict(risk_symbols=len(risk),mapped_UU=len(mappedUU),unresolved=len(set(risk)-mappedUU),nearest=sorted(nearest,key=lambda r:r['distance_bp'])[0]['UU_gene'],nearest_distance=sorted(nearest,key=lambda r:r['distance_bp'])[0]['distance_bp'],hits_within300kb=len(result['risk_hits_within300kb']))))
