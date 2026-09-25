from pathlib import Path
import csv,json,gzip,collections,urllib.request,xml.etree.ElementTree as ET,concurrent.futures,hashlib
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.23.0';O.mkdir(exist_ok=True)
rows=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'))
def fetch(r):
 url='https://www.ebi.ac.uk/ena/browser/api/xml/'+r['SRX'];b=urllib.request.urlopen(url,timeout=40).read();(O/(r['SRX']+'.xml')).write_bytes(b);e=ET.fromstring(b).find('EXPERIMENT');assert e.attrib['accession']==r['SRX']
 u='https://www.ebi.ac.uk/ena/portal/api/filereport?accession='+r['SRR']+'&result=read_run&fields=run_accession,submitted_format,submitted_ftp,submitted_md5,submitted_bytes&format=tsv';raw=urllib.request.urlopen(u,timeout=40).read();(O/(r['SRR']+'.submitted.tsv')).write_bytes(raw)
 return dict(sample=r['column'],experiment=r['SRX'],source=url,xml_sha256=hashlib.sha256(b).hexdigest(),protocol=e.findtext('.//LIBRARY_CONSTRUCTION_PROTOCOL'),layout=[x.tag for x in e.find('.//LIBRARY_LAYOUT')],title=e.findtext('TITLE'),submitted=list(csv.DictReader(raw.decode().splitlines(),delimiter='\t')))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:meta=list(pool.map(fetch,rows))
indexes=set('CTAGCT CGATGT TTAGGC TGACCA ACAGTG TATAAT CAGATC ACTTGA GATCAG TAGCTT GGCTAC CTTGTA'.split());out=[]
for r in rows:
 for mate in [1,2]:
  p=P/'06_v1.22.6/raw'/f"{r['SRR']}_{mate}.fastq.gz";n=0;bc=0;pt=0;joint=0;length=collections.Counter();phred=collections.Counter()
  with gzip.open(p,'rt') as f:
   for i in range(100000):
    h,s,plus,q=[f.readline().rstrip() for _ in range(4)];assert h.startswith('@') and plus.startswith('+') and len(s)==len(q)
    hit=s[:6] in indexes;t=s[6:26].count('T')>=16;n+=1;bc+=hit;pt+=t;joint+=hit and t;length[len(s)]+=1;phred.update(q)
  out.append(dict(sample=r['column'],mate=mate,n=n,known_inline_prefix_fraction=bc/n,polyT_80pct_positions7_26_fraction=pt/n,joint_index_polyT_fraction=joint/n,length_histogram=dict(length),phred_histogram={str(ord(k)-33):v for k,v in phred.items()}))
summary=dict(metadata=meta,read_structure=out,sampling='First 100000 reads in each of 18 raw ENA FASTQs; not random nor whole-file structure validation',manual='Takara authored 092618 manual pp8,21; mirror https://device.report/m/a1cff4f8333dfe571d210a1f84e734b6cdb0f611dc09ffa2d1e02c4b1109f9f6.pdf',limitations=['Metadata may be copied across repositories','Absence of inline structure could reflect preprocessing or a different library workflow; does not identify which','No raw-file modification or reinterpretation of barcode as UMI'])
(O/'protocol_structure_audit.json').write_text(json.dumps(summary,indent=2));print(json.dumps(dict(protocols=sorted(set(r['protocol'] for r in meta)),layouts=[r['layout'] for r in meta],R2_index_range=[min(x['known_inline_prefix_fraction'] for x in out if x['mate']==2),max(x['known_inline_prefix_fraction'] for x in out if x['mate']==2)],R2_polyT_range=[min(x['polyT_80pct_positions7_26_fraction'] for x in out if x['mate']==2),max(x['polyT_80pct_positions7_26_fraction'] for x in out if x['mate']==2)],R2_joint_range=[min(x['joint_index_polyT_fraction'] for x in out if x['mate']==2),max(x['joint_index_polyT_fraction'] for x in out if x['mate']==2)]),indent=2))
