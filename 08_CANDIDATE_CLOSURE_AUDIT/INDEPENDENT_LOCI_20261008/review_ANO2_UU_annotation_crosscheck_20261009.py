from pathlib import Path
import requests,gzip,hashlib,json,datetime,urllib.parse,collections
N=Path(__file__).resolve().parent;O=N/'ANO2_UU_annotation_crosscheck_20261009';O.mkdir(exist_ok=False)
base='https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9615/106/GCF_011100685.1_UU_Cfam_GSD_1.0/';name='GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gff.gz';s=requests.Session();resp=s.get(base+'md5checksums.txt',timeout=60);resp.raise_for_status();checks=resp.text;(O/'md5checksums.txt').write_text(checks);expected=[line.split()[0] for line in checks.splitlines() if line.split()[-1].lstrip('./')==name];assert len(expected)==1
p=O/name
with s.get(base+name,stream=True,timeout=120) as resp:
 resp.raise_for_status();headers=dict(resp.headers)
 with p.open('wb') as f:
  for block in resp.iter_content(1024*1024):f.write(block)
assert hashlib.md5(p.read_bytes()).hexdigest()==expected[0]
c='NC_049248.1';lo=6741835;hi=6742835;overlap=[];tss=[];counts=collections.Counter();directives=[];seq_region=None
with gzip.open(p,'rt') as f:
 for line in f:
  if line.startswith('#'):
   if line.startswith('##sequence-region '+c+' '):seq_region=line.strip()
   if line.startswith('#!'):directives.append(line.strip())
   continue
  x=line.rstrip().split('\t')
  if len(x)!=9 or x[0]!=c:continue
  kind=x[2];start=int(x[3])-1;end=int(x[4]);assert start>=0 and end>start;counts[kind]+=1
  attrs={k:urllib.parse.unquote(v) for item in x[8].split(';') if '=' in item for k,v in [item.split('=',1)]};record=dict(type=kind,start0=start,end0=end,strand=x[6],attributes=attrs)
  rna=kind=='mRNA' or kind.endswith('RNA') or kind in ('transcript','primary_transcript')
  if kind in ('gene','pseudogene','exon','CDS') or rna:
   if start<hi and end>lo:overlap.append(record)
  if rna and x[6] in ('+','-'):
   point=start if x[6]=='+' else end-1;d=lo-point if point<lo else (point-hi+1 if point>=hi else 0);tss.append(dict(record,tss0=point,minimum_distance_bp=d))
assert seq_region is not None and tss
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source=base+name,source_report='https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Canis_lupus_familiaris/106/',source_headers=headers,md5=expected[0],sha256=hashlib.sha256(p.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),UU_window=[c,lo,hi],sequence_region=seq_region,metadata=directives,chromosome_feature_counts=dict(counts),overlaps=overlap,nearest_mRNA_TSS=sorted([r for r in tss if r['type']=='mRNA'],key=lambda r:r['minimum_distance_bp'])[:20],nearest_all_RNA_TSS=sorted(tss,key=lambda r:r['minimum_distance_bp'])[:20],limits=['GFF transcript strand-aware annotated5prime boundaries, not measured alternative promoter activity.','Another annotation/assembly view of AR106, not independent experimental confirmation or functional isolation.','No risk-gene or miRNA clearance inferred from RNA-start distances alone; prior dedicated gates remain separate.'])
(O/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(state='completed',overlaps=len(overlap),nearest_mRNA_distance=result['nearest_mRNA_TSS'][0]['minimum_distance_bp'],nearest_all_RNA_distance=result['nearest_all_RNA_TSS'][0]['minimum_distance_bp'],nearest_mRNA=result['nearest_mRNA_TSS'][0]['attributes'].get('gene'))))
