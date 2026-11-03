from pathlib import Path
import json,hashlib
N=Path(__file__).resolve().parent;P=N.parent.parent
file=P/'05_SHIP/canine_all_genes_stranded.bed';genes=[]
for line in file.read_text().splitlines():
 f=line.split('\t');genes.append(dict(chrom=f[0],start0=int(f[1]),end0=int(f[2]),gene=f[3],strand=f[5]))
rows=[]
for w in json.loads((N/'next_queue.json').read_text()):
 if w['id'] not in ['bg1k_1059','bg1k_1266','bg1k_1900']:continue
 s=int(w['start0']);e=int(w['end0']);g=[r for r in genes if r['chrom']==w['chrom']];left=max((r for r in g if r['end0']<=s),key=lambda r:r['end0']);right=min((r for r in g if r['start0']>=e),key=lambda r:r['start0'])
 rows.append(dict(id=w['id'],chrom=w['chrom'],window_start0=s,window_end0=e,left=left,right=right,intergenic_start0=left['end0'],intergenic_end0=right['start0'],convergent=left['strand']=='+' and right['strand']=='-',interpretation='Flanking coding-gene interval is descriptive; not substituted for the original SHIP regional definition.'))
(N/'flanking_coding_context.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),results=rows),indent=2));print(json.dumps(rows,indent=2))
