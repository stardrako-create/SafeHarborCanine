from pathlib import Path
import csv,json
import numpy as np
import pyBigWig
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'06_v1.21.0'
rows=list(csv.DictReader((O/'raw_scoring.tsv').open(),delimiter='\t'))
b=pyBigWig.open(str(P/'06_v1.20.0/RRBS/methylation_weighted_mean.bw'))
c=pyBigWig.open(str(P/'06_v1.20.0/RRBS/cpg_coverage_frequency.bw'))
out=[]
for r in rows:
 if r['hard_veto']!='False':continue
 chrom,s,e=r['chrom'],int(r['start']),int(r['end'])
 v=b.values(chrom,s,e,numpy=True);n=c.values(chrom,s,e,numpy=True)
 mask=n>0
 out.append({'chrom':chrom,'start':s,'end':e,'genes':r['left_gene']+'/'+r['right_gene'],'reported_rrbs_mean':float(r['rrbs_mean']),'evidence_fraction':float(mask.mean()),'mean_on_evidence_bp':float(v[mask].mean()),'unobserved_bp_written_zero':int(((~mask)&(v==0)).sum())})
b.close();c.close()
with (O/'rrbs_missing_as_zero_audit.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
print(json.dumps(out,indent=2))
with (O/'AUDIT.md').open('a',encoding='utf-8') as f:
 f.write('\n## Achado adicional prioritário: média RRBS confunde ausência de dados com zero\n\nO builder escreve zero na média quando não há evidência. O scorer calcula a média sobre a janela inteira e apenas verifica se existe alguma cobertura algures na janela. A média reportada é por isso diluída pela fração de bases sem evidência. A estatística condicionada a coverage>0 foi medida diretamente nas quatro janelas e está em rrbs_missing_as_zero_audit.tsv. Esta é uma comparação diagnóstica, não uma recalibração do score: corrigir exige aplicar a mesma máscara aos candidatos e ao background. Todos os scores de low_methylation e rankings associados permanecem provisórios.\n')
 f.write('\n| locus | média atual (%) | fração com evidência | média nas bases com evidência (%) |\n|---|---:|---:|---:|\n')
 for r in out:f.write(f"| {r['genes']} | {r['reported_rrbs_mean']:.2f} | {r['evidence_fraction']:.3%} | {r['mean_on_evidence_bp']:.2f} |\n")
with (O/'CURRENT.md').open('a',encoding='utf-8') as f:f.write('\nCritical audit finding: RRBS mean includes unobserved bins as zero. Rankings remain provisional until coverage-aware candidate AND background scoring is implemented.\n')
