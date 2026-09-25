from pathlib import Path
import csv,json,gzip,hashlib,shutil,math
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.3';f=O/'GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz'
with gzip.open(f,'rt') as inp:reader=csv.DictReader(inp,delimiter='\t');columns=reader.fieldnames;rows=list(reader)
assert len(columns)==10 and len(rows)==14385 and len({r['Gene'] for r in rows})==14385
assert all(math.isfinite(float(r[c])) and float(r[c])>=0 for r in rows for c in columns[1:])
genes=['ANO2','NTF3','LOC119876429','LOC119872513','NPNT','TBCK','TSEN15','C7H1orf21','ARPC5','TET2','CD3D','CD3E','CD3G','CD247','LCK','CD4','CD8A','CD8B'];lookup={r['Gene']:r for r in rows};out=[]
for g in genes:out.append(dict(Gene=g,in_deposited_CPM_table=g in lookup,**{c:lookup.get(g,{}).get(c,'') for c in columns[1:]}))
with (O/'CAR_T_context_gene_inventory.tsv').open('w',newline='') as dest:w=csv.DictWriter(dest,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
source='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/GSE247355/suppl/GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz'
(O/'CAR_T_source_inventory.json').write_text(json.dumps(dict(accession='GSE247355',url=source,sha256=hashlib.sha256(f.read_bytes()).hexdigest(),gene_rows=len(rows),sample_columns=columns[1:],download_validated=True,limitations=['Deposited normalized CPM, not raw counts; do not feed into count-based DE testing','Sample prefixes B/E/M not yet independently resolved to donor metadata in this analysis','Missing genes in this processed table are not evidence of zero expression','RNA expression context does not measure accessibility or demonstrate safe integration']),indent=2))
note='''# Safe Harbor CAR-T — v1.22.3: profundidade do multiome e dados públicos CAR-T caninos

15 de setembro de 2026. **Há dados públicos diretamente de CAR-T caninas. A referência RNA GSE247355 foi agora identificada, descarregada e verificada. Além disso, a análise da profundidade mostra que os zeros do multiome tumoral têm pouca força para argumentar contra acessibilidade.** Nenhum destes resultados demonstra um safe harbor.

## Dados CAR-T caninos encontrados

[GSE247355](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247355), associado ao [estudo de Cao e colaboradores, 2024](https://pubmed.ncbi.nlm.nih.gov/38554158/), deposita RNA-seq de CAR-T B7-H3, CAR-T B7-H3/CXCR2 e células T controlo. Há nove amostras/colunas, organizadas em prefixos B, E e M. A identidade e condições dos dadores ainda não foram reconciliadas individualmente nesta análise; não se deve inferir automaticamente o desenho apenas dos nomes.

A matriz CPM pública foi descarregada: 14.385 linhas génicas e nove colunas numéricas, todas finitas e não negativas. Os reads brutos estão associados ao SRA. Foi criada uma tabela de presença e valores para genes próximos dos candidatos e marcadores T, como inventário inicial. Ausência de um gene nesta tabela processada não foi convertida em zero. Ainda não foi feita análise diferencial; CPM não deve ser entregue como contagens brutas a um modelo que exija counts.

Este dataset é diretamente útil para contexto de expressão em células CAR-T caninas. Não é ATAC-seq e não resolve diretamente se a janela de inserção está aberta. Também não demonstra neutralidade de integração nos nossos loci. A sua ausência nas iterações anteriores era uma lacuna da pesquisa, não ausência geral de dados CAR-T na literatura.

O estudo [Unedited allogeneic iNKT cells show extended persistence in MHC-mismatched canine recipients](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591065/) disponibiliza GSE229457 (scRNA-seq) e GSE229458 (Nanostring). Inclui investigação de iNKT e engenharia CAR; são recursos complementares de outra população celular. Não se deve classificar todo o material como CAR-T convencional, nem chamar RNA-seq ao painel Nanostring. As condições de cada amostra precisam de verificação antes de integração.

Para identidade celular canina independente, também foram identificados [GSE225599](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225599), com leucócitos de sete cães saudáveis e dez com osteossarcoma, e [GSE301630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301630), com seis cães saudáveis. São referências transcriptómicas, não ATAC CAR-T. Nesta pesquisa não foi confirmado um dataset público de ATAC-seq especificamente de CAR-T caninas. Isto descreve o resultado da pesquisa, não prova inexistência.

## Profundidade do grupo T no multiome usado

Recalculámos diretamente do H5 as métricas RNA/ATAC da matriz e confirmámos a concordância exata com o cache de todos os 5.969 barcodes. No conjunto QC de 5.849 núcleos, os 155 núcleos T consensuais representam 2,65% dos núcleos, mas apenas **1,07% das contagens da matriz de picos**. A mediana é 2.787 contagens por núcleo T, contra 6.064 nos restantes núcleos QC excluindo os 157 da união T.

Isto confirma uma diferença de profundidade da matriz. Não é total de fragmentos genómicos, FRiP, TSS enrichment ou um certificado de qualidade. Continuam por avaliar adequadamente nucleossomas, doublets e identidade independente.

## Os zeros são informativos?

Fizemos 2.000 amostragens descritivas, com semente 42, de grupos de 155 núcleos dos restantes QC, preservando a distribuição dos T em dez estratos de profundidade ATAC da matriz. Foram excluídos os 157 da união T do conjunto de comparação. Cada amostragem escolhe sem reposição dentro dos estratos; os grupos entre amostragens podem reutilizar núcleos. A profundidade agregada mediana destes grupos foi 559.883,5 contagens, próxima das 550.644 dos T.

Nas três janelas originais, os grupos de comparação com profundidade semelhante tiveram zero inícios de fragmentos em **63,45% (ANO2), 44,40% (LOC) e 41,45% (NPNT)** das amostragens. Portanto, zero nos T não é, por si só, um resultado surpreendente nesta escala de deteção. Isto reforça a conclusão já cautelosa: falta suporte positivo robusto, mas não há demonstração de cromatina fechada.

Nos controlos CD3D/CD3E/CD247, observaram-se 4/4/7 inícios nos T, enquanto as medianas dos grupos de comparação foram zero. É coerente com a anotação exploratória T, mas não valida pureza, identidade individual ou acessibilidade de todos os candidatos.

As distribuições de reamostragem são comparações condicionais dentro de um tumor, não intervalos de confiança entre dadores ou testes de significância biológica. O grupo de comparação inclui populações heterogéneas. A profundidade por picos é uma proxy e não corrige todos os vieses. A tabela cobre as 37 alternativas com mapeamento integral/contínuo e os sete intervalos originais/controlos; não estende a comparação integral às 27 alternativas parcialmente mapeadas.

## Prioridade atualizada

A próxima análise deve incorporar a expressão em CAR-T caninas de GSE247355 e usar referências caninas saudáveis para melhorar a identidade celular do multiome. O multiome tumoral permanece evidência auxiliar de acessibilidade esparsa. Nem os seus zeros nem a expressão de genes vizinhos devem virar um veto ou um passe automático.

Ficheiros: `depth_group_summary.tsv`, `depth_matched_local_counts.tsv`, `depth_review_summary.json`, matriz CPM GSE247355, `CAR_T_source_inventory.json`, `CAR_T_context_gene_inventory.tsv`. A pesquisa de fontes e o inventário não equivalem a uma análise completa dos novos estudos. A shortlist e o scorer global não foram alterados.
'''
(O/'CURRENT.md').write_text(note,encoding='utf-8');report=note+'\n\n---\n\n# Histórico anterior até v1.22.2\n\nAs fontes CAR-T e a avaliação da força dos zeros foram atualizadas acima.\n\n'+(P/'06_v1.22.2/RELATORIO_ATUALIZADO_2026-09-15.md').read_text(encoding='utf-8');(O/'RELATORIO_ATUALIZADO_2026-09-15.md').write_text(report,encoding='utf-8')
for name in ['depth_review_v1223.py','finish_depth_sources_v1223.py']:shutil.copy2(Path(name),O/name)
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.3 profundidade e fontes CAR-T - 2026-09-15';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');shared=notes/'Safe Harbor CAR-T - Relatorio atualizado - 2026-09-15.md';shutil.copy2(shared,O/'previous_Obsidian_report.md');shared.write_text(report,encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!important] Estado mais recente: v1.22.3\n> [[{name}]] — RNA-seq CAR-T canino GSE247355 descarregado; profundidade reduzida nos T do multiome limita a interpretação dos zeros. Sem alteração da shortlist.\n',1),encoding='utf-8')
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.22.3\n\n[Relatório](06_v1.22.3/RELATORIO_ATUALIZADO_2026-09-15.md).\n\nFonte RNA-seq CAR-T canina GSE247355 disponível e descarregada; análise completa pendente. Profundidade do multiome revista; zeros não demonstram cromatina fechada. Scoring global v1.21.8, sem safe harbor validado.\n',encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.3_Depth_and_Canine_CART_Sources','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'output_hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='output_hashes.json'},indent=2));print('v1.22.3 saved: depth review complete; canine CAR-T matrix downloaded and inventoried, full RNA analysis pending.')
