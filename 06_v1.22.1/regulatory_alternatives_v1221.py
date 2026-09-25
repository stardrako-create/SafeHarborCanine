from pathlib import Path
import csv,json,hashlib
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.1'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
windows=read(O/'review_windows.tsv');epic=read(P/'06_v1.21.5/epic_nonquiescent_ROS_blocks.tsv');sources={}
def bed(name):
 p=P/'05_SHIP'/name;sources[name]=hashlib.sha256(p.read_bytes()).hexdigest();return [v for line in p.read_text().splitlines() if line and not line.startswith('#') and len(v:=line.split('\t'))>=3]
genes=bed('canine_all_genes_stranded.bed');rna=bed('canine_lncRNA_smallRNA.bed');mirna=bed('canine_miRNA.bed');external=bed('ehsan_regulatory_elements_ROS.bed');risk={r['gene_symbol']:r['categories'] for r in read(P/'05_SHIP/canine_risk_genes.tsv')}
def distance(s,e,a,b):return max(s-b,a-e,0)
def union(xs):
 end=-1;total=0
 for s,e in sorted(xs):total+=max(0,e-max(s,end));end=max(end,e)
 return total
out=[];hits=[]
for r in windows:
 s=int(r['start']);e=int(r['end']);chrom=r['chrom'];local=[g for g in genes if g[0]==chrom];body=min(local,key=lambda g:distance(s,e,int(g[1]),int(g[2])));tss=min(local,key=lambda g:distance(s,e,int(g[1]) if g[5]=='+' else int(g[2])-1,(int(g[1]) if g[5]=='+' else int(g[2])-1)+1));pos=int(tss[1]) if tss[5]=='+' else int(tss[2])-1
 risklocal=[g for g in local if g[3] in risk];rg=min(risklocal,key=lambda g:distance(s,e,int(g[1]),int(g[2]))) if risklocal else None
 states=[a for a in epic if a['genes']==r['genes'] and int(a['ros_start'])<e and int(a['ros_end'])>s];reg=[a for a in states if 1<=int(a['state'])<=7]
 regbp=union([(max(s,int(a['ros_start'])),min(e,int(a['ros_end']))) for a in reg]);assert regbp==int(r['epic_promoter_enhancer_bp_union'])==0
 result=dict(window_id=r['window_id'],nearest_gene_body=body[3],gene_body_gap_bp=distance(s,e,int(body[1]),int(body[2])),nearest_gene_5prime=tss[3],gene_5prime_gap_bp=distance(s,e,pos,pos+1),nearest_risk_gene=rg[3] if rg else '',risk_gene_categories=risk.get(rg[3],'') if rg else '',risk_gene_body_gap_bp=distance(s,e,int(rg[1]),int(rg[2])) if rg else '',epic_promoter_enhancer_bp=regbp,epic_observed_nonquiescent_states=';'.join(str(x) for x in sorted({int(a['state']) for a in states})))
 for name,data in [('lncrna_smallrna',rna),('mirna',mirna),('external_regulatory',external)]:
  intervals=[(max(s,int(v[1])),min(e,int(v[2]))) for v in data if v[0]==chrom and int(v[1])<e and int(v[2])>s];result[name+'_overlap_bp']=union(intervals)
 for a in states:hits.append(dict(window_id=r['window_id'],tissue=a['tissue'],state=a['state'],overlap_bp=min(e,int(a['ros_end']))-max(s,int(a['ros_start']))))
 out.append(result)
for name,rows in [('regulatory_context.tsv',out),('epic_state_hits.tsv',hits)]:
 if rows:
  with (O/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(O/'regulatory_provenance.json').write_text(json.dumps(dict(source_sha256=sources,epic_blocks_sha256=hashlib.sha256((P/'06_v1.21.5/epic_nonquiescent_ROS_blocks.tsv').read_bytes()).hexdigest(),limitations=['Distances use existing gene annotation and one gene 5-prime boundary, not exhaustive alternative transcript TSS','Gap convention: touching half-open intervals have zero gap','No new threshold chosen; descriptive context only','EpiC mapped eleven tissues, no purified T; no overlap does not establish regulatory inactivity','Risk gene list is an existing project list, not an exhaustive oncogenicity assessment']),indent=2));print(json.dumps(out,indent=2))
