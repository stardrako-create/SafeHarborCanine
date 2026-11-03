from pathlib import Path
import csv,json,hashlib,collections
O=Path(__file__).parent;P=O.parent.parent
windows=json.loads((O/'UU_projection/inputs.json').read_text())['windows']
names=['canine_all_genes.bed','canine_miRNA.bed','canine_lncRNA_smallRNA.bed','ehsan_regulatory_elements_ROS.bed']
beds={}
for name in names:
 b=collections.defaultdict(list)
 for l in (P/'05_SHIP'/name).read_text().splitlines():
  if not l or l.startswith('#'):continue
  v=l.split();b[v[0]].append((int(v[1]),int(v[2]),v[3:] ))
 beds[name]=b
risk={r['gene_symbol'] for r in csv.DictReader((P/'05_SHIP/canine_risk_genes.tsv').open(),delimiter='\t')}
out=[]
for w in windows:
 c=w['chrom'];s=int(w['start']);e=int(w['end']);r={ 'id':f'{c}_{s}_{e}'}
 for name in names:
  hits=[(a,b,x) for a,b,x in beds[name].get(c,[]) if a<e and b>s]
  r[name+'_overlaps']=hits
 genes=beds['canine_all_genes.bed'].get(c,[])
 targets=[(max(a-e,s-b,0),x[0]) for a,b,x in genes if x and x[0] in risk]
 r['nearest_catalogue_risk_gene']=min(targets) if targets else None
 mir=[max(a-e,s-b,0) for a,b,x in beds['canine_miRNA.bed'].get(c,[])]
 r['nearest_miRNA_gap']=min(mir) if mir else None
 out.append(r)
(O/'direct_annotation_review_20261008.json').write_text(json.dumps(dict(results=out,source_hashes={n:hashlib.sha256((P/'05_SHIP'/n).read_bytes()).hexdigest() for n in names+['canine_risk_genes.tsv']},limitations=['Catalogue annotations only; not absence of unannotated regulation']),indent=2))
print(json.dumps(next(x for x in out if x['id']=='NC_051831.1_39708074_39709074')))