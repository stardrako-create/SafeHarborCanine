from pathlib import Path
import json,csv,requests,hashlib,datetime
N=Path(__file__).resolve().parent;O=N/'ANO2_local_ortholog_risk_review_20261009';O.mkdir(exist_ok=False);inventory=json.loads((N/'ANO2_local_gene_inventory_20261009.json').read_text());genes={g['GeneIDs'][0]:g for r in inventory['results'] for g in r['annotated_genes_within300kb']};assert len(genes)==6
H=N/'risk_HGNC_rename_audit_20261009';hgnc={r['hgnc_id']:r for r in csv.DictReader((H/'hgnc_complete_set.txt').open(),delimiter='\t')};links=json.loads((H/'risk_symbol_links.json').read_text());riskids={}
for r in links:
 for key in r['HGNC_ids']:
  gid=hgnc[key]['entrez_id']
  if gid:riskids.setdefault(gid,[]).append(r['original_risk_symbol'])
session=requests.Session();results=[];status=dict(state='running',completed=[])
def save():(O/'status.json').write_text(json.dumps(status,indent=2))
save()
for gid,g in genes.items():
 url=f'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gid}/orthologs';record=dict(canine_gene=g,source_url=url)
 if gid=='486731':
  source=N/'local_gene_ortholog_probe_20261009/raw_response.txt';raw=source.read_bytes();code=200;record['reused_probe_response_sha256']=hashlib.sha256(raw).hexdigest()
 else:
  r=session.get(url,params={'page_size':1000},timeout=60);raw=r.content;code=r.status_code;record['request_url']=r.url;record['response_headers']=dict(r.headers)
 (O/(gid+'.response.json')).write_bytes(raw);record['raw_sha256']=hashlib.sha256(raw).hexdigest();record['http_status']=code
 if code==200:
  data=json.loads(raw);rr=data.get('reports',[]);record.update(total_count=data.get('total_count'),returned_count=len(rr),next_page_token=data.get('next_page_token'))
  assert not data.get('next_page_token'),'Ortholog page truncated; no negative conclusion allowed'
  assert int(data.get('total_count',len(rr)))==len(rr),'Ortholog count mismatch'
  humans=[x['gene'] for x in rr if str(x['gene']['tax_id'])=='9606'];record['human_orthologs']=[dict(gene_id=h['gene_id'],symbol=h['symbol'],type=h.get('type'),gene_groups=h.get('gene_groups'),matching_frozen_risk_symbols=riskids.get(str(h['gene_id']),[])) for h in humans];record['status']='human_ortholog_compared_by_GeneID' if humans else 'no_human_ortholog_in_returned_dataset_not_risk_absence'
 else:record['status']='ortholog_metadata_unavailable_not_risk_absence';record['error_response']=raw.decode(errors='replace')[:1000]
 results.append(record);status['completed'].append(dict(gene=g['gene'],status=record['status']));save()
receipt=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),inventory_sha256=hashlib.sha256((N/'ANO2_local_gene_inventory_20261009.json').read_bytes()).hexdigest(),frozen_risk_HGNC_links_sha256=hashlib.sha256((H/'risk_symbol_links.json').read_bytes()).hexdigest(),risk_entrez_geneids=len(riskids),results=results,known_human_risk_GeneID_matches=sum(bool(h['matching_frozen_risk_symbols']) for r in results for h in r.get('human_orthologs',[])),limits=['Only annotated local genes within300kb, using NCBI ortholog metadata and frozen risk-list HGNC GeneIDs; no gene-role expansion.','Missing human ortholog for local lncRNA is unknown risk-role correspondence, not absence of biological role or risk.','NCBI computational orthology is not independent biological confirmation or proof of functional safety.','Source taxon filtering in the initial NTF3 probe was not honored; entire returned set was count-checked and human records filtered locally.'])
(O/'receipt.json').write_text(json.dumps(receipt,indent=2));status.update(state='completed',known_human_risk_GeneID_matches=receipt['known_human_risk_GeneID_matches']);save();print(json.dumps(status))
