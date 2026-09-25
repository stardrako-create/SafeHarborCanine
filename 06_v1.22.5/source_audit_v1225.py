from pathlib import Path
import csv,gzip,json,re,hashlib,xml.etree.ElementTree as ET,urllib.request,concurrent.futures
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.5'
with gzip.open(O/'Canis_familiaris.gene_info.gz','rt') as f:info=list(csv.DictReader(f,delimiter='\t'))
assert all(r['#tax_id']=='9615' for r in info)
with gzip.open(P/'06_v1.22.3/GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz','rt') as f:matrix=list(csv.DictReader(f,delimiter='\t'))
symbols={r['Gene'] for r in matrix};selected=sorted({r['gene'] for r in csv.DictReader((P/'06_v1.22.4/CAR_T_gene_expression_context.tsv').open(),delimiter='\t')});aliases={}
for r in info:
 for name in {r['Symbol'],r['Symbol_from_nomenclature_authority'],*r['Synonyms'].split('|'),'LOC'+r['GeneID']} - {'-',''}:aliases.setdefault(name,[]).append(r)
out=[]
for g in selected:
 matches={r['GeneID']:r for r in aliases.get(g,[])};candidate={}
 for gid,r in matches.items():
  names={r['Symbol'],r['Symbol_from_nomenclature_authority'],*r['Synonyms'].split('|'),'LOC'+gid}-{'-',''}
  for n in names&symbols:candidate[n]=sorted({x['GeneID'] for x in aliases.get(n,[])})
 unambiguous=len(matches)==1 and len(candidate)==1 and next(iter(candidate.values()),[])==list(matches)
 out.append(dict(requested_symbol=g,exact_in_CPM=g in symbols,NCBI_gene_ids=';'.join(matches),current_NCBI_symbols=';'.join(r['Symbol'] for r in matches.values()),CPM_alias_candidates=';'.join(sorted(candidate)),alias_symbol_to_gene_ids=json.dumps(candidate),unambiguous_alias_recovery=bool(g not in symbols and unambiguous),status='exact_present' if g in symbols else 'recoverable_by_current_unique_gene_alias' if unambiguous else 'no_CPM_alias_match' if not candidate else 'ambiguous_requires_review'))
with (O/'gene_alias_audit.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
root=ET.parse(O/'PMC10981605.xml').getroot();supp=[]
for e in root.iter('supplementary-material'):
 supp.append(dict(id=e.attrib.get('id'),links=[v for child in e.iter() for k,v in child.attrib.items() if k.endswith('href')]))
paragraphs=[' '.join(''.join(e.itertext()).split()) for e in root.iter() if e.tag in ['p','caption']]
relevant=[t for t in paragraphs if any(x in t.lower() for x in ['3 individual dogs','3 different dogs','htseq','ruvseq','day 10','10 days of expansion'])]
(O/'paper_method_evidence.json').write_text(json.dumps(dict(paper='https://doi.org/10.1007/s00262-024-03642-4',fulltext_source='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10981605/fullTextXML',relevant_passages=relevant,supplementary_links=supp),indent=2))
meta=list(csv.DictReader((P/'06_v1.22.4/sample_metadata_audit.tsv').open(),delimiter='\t'))
def fetch(r):
 srx=re.search(r'SRX\d+',r['relations']).group();url='https://www.ebi.ac.uk/ena/portal/api/filereport?accession='+srx+'&result=read_run&fields=run_accession,experiment_accession,sample_accession,sample_alias,experiment_alias,library_name,fastq_ftp,fastq_bytes&format=tsv';file=O/(srx+'.ENA.tsv')
 try:
  with urllib.request.urlopen(url,timeout=30) as f:data=f.read()
  file.write_bytes(data);records=list(csv.DictReader(data.decode().splitlines(),delimiter='\t'))
  return dict(GSM=r['accession'],column=r['column'],SRX=srx,url=url,records=records,error=None)
 except Exception as e:return dict(GSM=r['accession'],column=r['column'],SRX=srx,url=url,records=[],error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:runs=list(pool.map(fetch,meta))
(O/'ENA_raw_read_inventory.json').write_text(json.dumps(runs,indent=2))
summary=dict(selected_genes=len(out),missing_exact=sum(not r['exact_in_CPM'] for r in out),recovered_by_unique_current_alias=[r for r in out if r['unambiguous_alias_recovery']],unresolved_missing=[r['requested_symbol'] for r in out if not r['exact_in_CPM'] and not r['unambiguous_alias_recovery']],ENA_runs=sum(len(r['records']) for r in runs),ENA_errors=[r['error'] for r in runs if r['error']],fastq_bytes=sum(int(v) for r in runs for x in r['records'] for v in x.get('fastq_bytes','').split(';') if v.isdigit()),limitations=['Current NCBI gene_info is not necessarily the annotation release used by authors','Absence after alias search does not determine filtering versus unannotated versus low expression','Raw-read sample aliases are additional metadata, not genotype-based identity verification','No raw FASTQ downloaded or reprocessed in this audit'])
(O/'source_audit_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(json.dumps(supp))
