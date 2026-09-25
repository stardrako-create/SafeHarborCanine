from pathlib import Path
import json,hashlib,shutil,csv
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.8';s=json.loads((O/'concordance_summary.json').read_text())
assert len(s['assessments'])==6 and all(x['identity_assignment'] and x['diagonal_best_count']==9 for x in s['assessments'])
table='| Genes comparados | Métrica | Correlação entre rótulos iguais | Melhor correspondência individual |\n|---|---|---:|---:|\n'
for a in s['assessments']:table+=f"| {a['genes']} ({a['gene_set']}) | {a['metric']} | {a['diagonal_min']:.4f}–{a['diagonal_max']:.4f} | 9/9 |\n"
note=f'''# v1.22.8 — concordância da reanálise com o RNA publicado

19 setembro 2026. **As nove amostras reprocessadas correspondem melhor às colunas publicadas com o mesmo rótulo em todas as seis avaliações.** Não foi encontrada evidência de troca de colunas entre a nossa reanálise e a matriz publicada. Isto não verifica identidade física dos cães nem exclui erros de origem partilhados pelos dois conjuntos.

## Dados e método

Comparada a matriz não orientada corrigida v1.22.7 com GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz. Dos 14.385 símbolos publicados, 11.976 têm correspondência exata e única no campo gene_name do índice ROS. Os 2.409 restantes foram excluídos desta comparação, com razões guardadas; isto não significa genes ausentes no animal. Não foram forçadas equivalências por aliases nem afirmada identidade das definições exónicas entre assemblies.

CPM reprocessado usa como denominador a soma de pares atribuídos a genes. Compararam-se todos os 81 pares de colunas: Spearman da abundância e Pearson de log2(CPM+1) após subtrair, em cada gene e separadamente em cada matriz, a média das nove amostras. A segunda análise reduz o efeito das diferenças de abundância entre genes e avalia melhor a correspondência dos perfis entre amostras. O +1 serve exclusivamente para esta transformação descritiva, sem teste diferencial.

Repetiu-se em genes com CPM≥1 em pelo menos três amostras de cada matriz (9.690 genes) e CPM≥5 (8.357 genes). São análises de sensibilidade da concordância, não novos filtros de candidatos. A atribuição global um-para-um que maximiza a soma de correlações também preserva todos os rótulos nas seis avaliações.

{table}

## O que fica esclarecido

A ligação operacional entre runs, nomes de amostras da reanálise e colunas publicadas tem suporte empírico. Os resultados são compatíveis com diferenças esperadas de anotação e método de contagem; a correlação elevada não prova equivalência quantitativa gene a gene. Esta comparação usa os mesmos dados biológicos, portanto não é replicação independente.

## O que permanece aberto

Os quatro conflitos source/tissue do GEO continuam preservados em metadata_conflicts_preserved.tsv. Não se alteraram rótulos nem se declarou resolvida a identidade física. A confirmação por genótipos ou esclarecimento do depositante seria evidência adicional diferente desta concordância; não foi enviado contacto.

O protocolo declarado SMART-Seq v4 3’ DE, a história de processamento dos FASTQs e a orientação continuam a exigir reconciliação. A consistência encontrada não justifica inventar trimming de índices/UMIs nem transformar o equilíbrio forward/reverse em confirmação definitiva do kit.

Antes de inferência: documentar normalização, correspondência génica e desenho pareado de três cães. As nove bibliotecas não são nove animais independentes. Não foram calculados p-values nem declarada ausência de efeito pela simples proximidade de expressão.

Para safe harbor: a evidência RNA fornece contexto em CAR-T; continua a faltar demonstração de acessibilidade e desempenho/segurança da integração nos intervalos candidatos. Nenhum locus foi promovido ou excluído apenas por esta análise.

## Artefactos

concordance_summary.json; cross_sample_concordance.tsv (486 correlações); exact_symbol_mapping.tsv; unmapped_published_symbols.tsv; metadata_conflicts_preserved.tsv. Scripts e hashes preservados nesta pasta. Última recontagem validada: v1.22.7. Não foi iniciado novo download nem alinhamento; acompanhamento da recontagem permanece pausado porque terminou.
'''
(O/'RELATORIO_CONCORDANCIA_2026-09-19.md').write_text(note,encoding='utf-8');(O/'CURRENT.md').write_text(note,encoding='utf-8')
for n in ['concordance_v1228.py','finish_concordance_v1228.py']:shutil.copy2(n,O/n)
paths=[P/'06_v1.22.7/gene_counts_unstranded.tsv',P/'06_v1.22.7/STAR_index/geneInfo.tab',P/'06_v1.22.3/GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz',P/'06_v1.22.4/sample_metadata_audit.tsv']
(O/'input_sha256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2))
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.22.8 — concordância RNA revista\n\nRecontagem v1.22.7 concluída e validada. Concordância v1.22.8: 9/9 perfis correspondem aos mesmos rótulos publicados em seis avaliações. Conflitos de metadados de origem e validação biológica continuam em aberto.\n\n[Relatório de concordância](06_v1.22.8/RELATORIO_CONCORDANCIA_2026-09-19.md) · [Recontagem corrigida](06_v1.22.7/RELATORIO_RECONTAGEM_CORRIGIDA_2026-09-19.md).\n',encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.8 concordancia RNA publicado - 2026-09-19';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!info] v1.22.8 — concordância com RNA publicado\n> [[{name}]]. 9/9 rótulos concordantes; conflitos de origem não resolvidos por esta comparação.\n',1),encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.8_PUBLISHED_RNA_CONCORDANCE','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print('Concordance report saved in project, Obsidian and MemPalace.')
