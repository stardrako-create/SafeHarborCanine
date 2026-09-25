from pathlib import Path
import csv,json,collections,hashlib,datetime
P=Path(__file__).parent.parent; O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922'
def bed(p):
 out=collections.defaultdict(list)
 for l in p.open():
  if l.startswith('#') or not l.strip():continue
  x=l.rstrip().split('\t');out[x[0]].append((int(x[1]),int(x[2]),x[3:]))
 return out
risk={r['gene_symbol'] for r in csv.DictReader((P/'05_SHIP/canine_risk_genes.tsv').open(),delimiter='\t')};genes=bed(P/'05_SHIP/canine_all_genes.bed');tads=bed(P/'05_SHIP/canine_tad_intervals.bed')
rows=list(csv.DictReader((O/'transcript_and_boundary_review.tsv').open(),delimiter='\t'));out=[]
for r in rows:
 c=r['chrom'];s=int(r['start0']);e=int(r['end0']);regions=[(a,b) for a,b,_ in tads.get(c,[]) if a<=s and b>=e]
 hits=sorted({extra[0] for a,b in regions for g,h,extra in genes.get(c,[]) if extra and extra[0] in risk and g<b and h>a})
 status='no_containing_TAD_interval_unresolved' if not regions else 'risk_catalogue_overlap' if hits else 'no_risk_catalogue_hit_in_containing_intervals'
 out.append(dict(id=r['id'],group=r['group'],chrom=c,start0=s,end0=e,containing_intervals=json.dumps(regions),risk_genes=';'.join(hits),status=status,scope='Existing interval proxy only; no T-cell contact validation or multi-scale audit'))
with (O/'TAD_risk_context_review.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
summary=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reviewed=len(out),counts=dict(collections.Counter(r['status'] for r in out)),sources_sha256={n:hashlib.sha256((P/'05_SHIP'/n).read_bytes()).hexdigest() for n in ['canine_tad_intervals.bed','canine_all_genes.bed','canine_risk_genes.tsv']},limitations=['Containing interval proxy only; no absence-of-risk conclusion for uncovered windows','No experimental T-cell TAD evidence','Multi-scale sensitivity remains pending for these control windows'])
(O/'TAD_risk_context_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
