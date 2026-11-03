from pathlib import Path
import requests,gzip,csv,json,collections,urllib.parse,hashlib,datetime,io,traceback
N=Path(__file__).resolve().parent;P=N.parent.parent;O=N/'risk_NCBI_canine_GeneID_followup_20261009';O.mkdir(exist_ok=False);url='https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Canis_lupus_familiaris.gene_info.gz';status=dict(state='running',stage='download_current_canine_gene_info');(O/'status.json').write_text(json.dumps(status))
try:
 resp=requests.get(url,timeout=120);resp.raise_for_status();raw=resp.content;(O/'Canis_lupus_familiaris.gene_info.gz').write_bytes(raw);data=gzip.decompress(raw).decode();rows=list(csv.DictReader(data.splitlines(),delimiter='\t'));assert rows and '#tax_id' in rows[0] and all(r['#tax_id']=='9615' for r in rows)
 H=N/'risk_HGNC_rename_audit_20261009';previous=json.loads((H/'receipt.json').read_text());remaining=set(previous['unmatched_union']);links=json.loads((H/'risk_symbol_links.json').read_text());look=collections.defaultdict(set)
 for r in links:
  if r['original_risk_symbol'] not in remaining:continue
  for spelling in set(r['accepted_same_record_symbol_spellings'])|{r['original_risk_symbol']}:look[spelling].add(r['original_risk_symbol'])
 ids=collections.defaultdict(set);info=[]
 for r in rows:
  names={r['Symbol'],r.get('Symbol_from_nomenclature_authority','')}|set(r['Synonyms'].split('|'));matched=set().union(*(look.get(name,set()) for name in names))
  if not matched:continue
  ids[r['GeneID']]|=matched;info.append(dict(matched_risk_symbols=sorted(matched),matched_spellings=sorted(names&set(look)),gene_info_record=r))
 (O/'canine_GeneID_candidates.json').write_text(json.dumps(info,indent=2));sources=[P/'01_referencia/ROS_Cfam_1.0/genomic.gff',N/'ANO2_UU_annotation_crosscheck_20261009/GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gff.gz'];results=[];placed=set();seenids=set()
 for path,c,lo,hi in [(sources[0],'NC_051831.1',39708074,39709074),(sources[1],'NC_049248.1',6741835,6742835)]:
  mappings=[];near=[];opener=gzip.open if path.suffix=='.gz' else open
  with opener(path,'rt') as f:
   for line in f:
    if line.startswith('#'):continue
    x=line.rstrip().split('\t')
    if len(x)!=9 or x[2]!='gene':continue
    a={k:urllib.parse.unquote(v) for item in x[8].split(';') if '=' in item for k,v in [item.split('=',1)]};gids={item.split(':',1)[1] for item in a.get('Dbxref','').split(',') if item.startswith('GeneID:')};matched=set().union(*(ids.get(gid,set()) for gid in gids))
    if not matched:continue
    seenids|=gids&set(ids);placed|=matched;start=int(x[3])-1;end=int(x[4]);record=dict(annotation_gene=a.get('gene'),GeneIDs=sorted(gids),matched_risk_symbols=sorted(matched),chrom=x[0],start0=start,end0=end,attributes=a);mappings.append(record)
    if x[0]==c:
     d=lo-end+1 if end<=lo else (start-hi+1 if start>=hi else 0);near.append(dict(record,distance_bp=d))
  results.append(dict(source=str(path),mapped_gene_records=mappings,nearest_risk_genes=sorted(near,key=lambda r:r['distance_bp'])[:20],hits_within300kb=[r for r in near if r['distance_bp']<300000]))
 found=set().union(*ids.values()) if ids else set();result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_url=url,source_headers=dict(resp.headers),gene_info_sha256=hashlib.sha256(raw).hexdigest(),source_HGNC_receipt_sha256=hashlib.sha256((H/'receipt.json').read_bytes()).hexdigest(),unresolved_input_symbols=len(remaining),symbols_with_current_canine_GeneID=len(found),new_symbols_placed_in_frozen_annotations=len(placed),newly_placed_symbols=sorted(placed),still_unplaced_symbols=sorted(remaining-placed),current_canine_GeneIDs_absent_from_frozen_annotation=sorted(set(ids)-seenids),assemblies=results,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),limits=['Current NCBI canine taxon9615 symbol/synonym records linked by exact GeneID to frozen AR106 coordinates.','Human risk-role transfer remains the existing proxy, not independently validated canine risk function.','Unplaced or absent GeneIDs do not prove gene absence; no coordinate or orthologue invented.','One-base differences from earlier gap reports reflect nearest-nucleotide vs interinterval-gap conventions, with unchanged300kb criterion and no borderline candidate cases.'])
 (O/'receipt.json').write_text(json.dumps(result,indent=2));status.update(state='completed',stage='coordinate_review_complete',newly_placed=len(placed),still_unplaced=len(remaining-placed),hits_within300kb=sum(len(r['hits_within300kb']) for r in results));print(json.dumps(status))
except Exception as e:status.update(state='failed',error=str(e));(O/'error.txt').write_text(traceback.format_exc());print(json.dumps(status))
finally:(O/'status.json').write_text(json.dumps(status,indent=2))
