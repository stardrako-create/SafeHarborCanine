from pathlib import Path
import requests,csv,json,gzip,urllib.parse,collections,hashlib,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent;O=N/'risk_HGNC_rename_audit_20261009';O.mkdir(exist_ok=False);url='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt';resp=requests.get(url,timeout=120);resp.raise_for_status();(O/'hgnc_complete_set.txt').write_bytes(resp.content);rows=list(csv.DictReader(resp.text.splitlines(),delimiter='\t'));assert rows and {'hgnc_id','symbol','status','prev_symbol','alias_symbol'}<=set(rows[0]);hgnc={r['hgnc_id']:r for r in rows if r['status']=='Approved'};idx={field:collections.defaultdict(set) for field in ('symbol','prev_symbol','alias_symbol')}
for key,r in hgnc.items():
 for field in idx:
  for name in r[field].split('|'):
   if name:idx[field][name].add(key)
def resolve(name):
 for field in ('symbol','prev_symbol','alias_symbol'):
  if name in idx[field]:return field,idx[field][name]
 return None,set()
riskfile=P/'05_SHIP/canine_risk_genes.tsv';risks={r['gene_symbol']:r['categories'] for r in csv.DictReader(riskfile.open(),delimiter='\t')};lookup=collections.defaultdict(set);links=[]
for original,categories in risks.items():
 field,ids=resolve(original);safe=set();unsafe=[]
 for key in ids:
  r=hgnc[key]
  for symbol in {r['symbol']}|set(r['prev_symbol'].split('|'))|set(r['alias_symbol'].split('|')):
   if not symbol:continue
   _,other=resolve(symbol)
   if other=={key}:safe.add(symbol)
   else:unsafe.append(dict(symbol=symbol,HGNC_ids=sorted(other)))
 for symbol in safe:lookup[symbol].add(original)
 lookup[original].add(original)
 links.append(dict(original_risk_symbol=original,categories=categories,HGNC_match_field=field,HGNC_ids=sorted(ids),ambiguous_original=len(ids)>1,accepted_same_record_symbol_spellings=sorted(safe),ambiguous_spellings_not_used=unsafe))
(O/'risk_symbol_links.json').write_text(json.dumps(links,indent=2));sources=[P/'01_referencia/ROS_Cfam_1.0/genomic.gff',N/'ANO2_UU_annotation_crosscheck_20261009/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gff.gz'];assemblies=[];allmatched=set()
for path,chrom,lo,hi in [(sources[0],'NC_051831.1',39708074,39709074),(sources[1],'NC_049248.1',6741835,6742835)]:
 matched=set();near=[];opener=gzip.open if path.suffix=='.gz' else open
 with opener(path,'rt') as f:
  for line in f:
   if line.startswith('#'):continue
   x=line.rstrip().split('\t')
   if len(x)!=9 or x[2]!='gene':continue
   a={k:urllib.parse.unquote(v) for item in x[8].split(';') if '=' in item for k,v in [item.split('=',1)]};names={a.get('gene',''),a.get('Name','')}|set(a.get('gene_synonym','').split(','));hit=set().union(*(lookup.get(name,set()) for name in names));matched|=hit
   if x[0]==chrom and hit:
    start=int(x[3])-1;end=int(x[4]);d=lo-end+1 if end<=lo else (start-hi+1 if start>=hi else 0);near.append(dict(gene=a.get('gene'),start0=start,end0=end,matched_original_risk_symbols=sorted(hit),distance_bp=d,attributes=a))
 allmatched|=matched;assemblies.append(dict(source=str(path),window=[chrom,lo,hi],matched_risk_symbols=len(matched),unmatched_symbols=sorted(set(risks)-matched),nearest=sorted(near,key=lambda r:r['distance_bp'])[:20],hits_within300kb=[r for r in near if r['distance_bp']<300000]))
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),HGNC_url=url,HGNC_source_headers=dict(resp.headers),HGNC_sha256=hashlib.sha256(resp.content).hexdigest(),official_field_definitions='https://www.genenames.org/help/statistics-and-downloads/',risk_list_symbols=len(risks),matched_union=len(allmatched),unmatched_union=sorted(set(risks)-allmatched),HGNC_unknown_risk_symbols=[r['original_risk_symbol'] for r in links if not r['HGNC_ids']],HGNC_ambiguous_risk_symbols=[r for r in links if r['ambiguous_original']],assemblies=assemblies,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [riskfile]+sources},script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limits=['Original risk membership/categories unchanged; a spelling/rename sensitivity using official human HGNC records, not a new risk list.','Canonical human symbol matching remains a canine orthology proxy; no orthology is proved solely by a human alias.','Ambiguous accepted original symbols are exposed and not a unique orthology assertion; ambiguous additional spellings excluded.','Unmatched symbols and missing annotation not classified as absence of risk; not biological confirmation.'])
(O/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(state=result['state'],risk_symbols=len(risks),matched_union=len(allmatched),unmatched=len(result['unmatched_union']),ambiguous_originals=len(result['HGNC_ambiguous_risk_symbols']),assemblies=[dict(source=r['source'],matched=r['matched_risk_symbols'],within300kb=len(r['hits_within300kb']),nearest=r['nearest'][0]['gene'],distance=r['nearest'][0]['distance_bp']) for r in assemblies])))
