from pathlib import Path
import gzip,json,hashlib,csv,urllib.parse,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent;S=N/'ANO2_UU_annotation_crosscheck_20261009';gff=S/'GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gff.gz';riskfile=P/'05_SHIP/canine_risk_genes.tsv';risk={r['gene_symbol']:r['categories'] for r in csv.DictReader(riskfile.open(),delimiter='\t')};matched=set();rgenes=[];mirna=[];lo=6741835;hi=6742835
with gzip.open(gff,'rt') as f:
 for line in f:
  if line.startswith('#'):continue
  x=line.rstrip().split('\t')
  if len(x)!=9 or x[2] not in ('gene','miRNA'):continue
  a={k:urllib.parse.unquote(v) for item in x[8].split(';') if '=' in item for k,v in [item.split('=',1)]};name=a.get('gene');start=int(x[3])-1;end=int(x[4]);d=lo-end+1 if end<=lo else (start-hi+1 if start>=hi else 0)
  if x[2]=='gene' and name in risk:
   matched.add(name)
   if x[0]=='NC_049248.1':rgenes.append(dict(gene=name,categories=risk[name],start0=start,end0=end,strand=x[6],distance_bp=d,attributes=a))
  if x[0]=='NC_049248.1' and x[2]=='miRNA':mirna.append(dict(start0=start,end0=end,strand=x[6],distance_bp=d,attributes=a))
assert rgenes and mirna
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),UU_window=['NC_049248.1',lo,hi],sources_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [gff,riskfile]},risk_symbols=len(risk),matched_risk_symbols_genomewide=len(matched),unmatched_risk_symbols=sorted(set(risk)-matched),nearest_risk_genes=sorted(rgenes,key=lambda r:r['distance_bp'])[:20],nearest_miRNAs=sorted(mirna,key=lambda r:r['distance_bp'])[:20],risk_genes_within_300kb=[r for r in rgenes if r['distance_bp']<300000],miRNAs_within_300kb=[r for r in mirna if r['distance_bp']<300000],limits=['Gene-symbol membership using the frozen project risk list, no claim all aliases/genes were mapped; unmatched symbols explicitly retained.','Distances use annotated extents, not functional effects or experimentally measured regulation.','UU annotation crosscheck supplements ROS gates; not proof of absence of risk outside annotation.'])
p=N/'ANO2_UU_risk_miRNA_crosscheck_20261009.json';assert not p.exists();p.write_text(json.dumps(result,indent=2));print(json.dumps(dict(risk_symbols=len(risk),matched=len(matched),nearest_risk=result['nearest_risk_genes'][0]['gene'],risk_distance=result['nearest_risk_genes'][0]['distance_bp'],miRNA_distance=result['nearest_miRNAs'][0]['distance_bp'],risk_under300kb=len(result['risk_genes_within_300kb']),miRNA_under300kb=len(result['miRNAs_within_300kb']))))
