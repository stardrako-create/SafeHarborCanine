# v1.22.8 — concordância da reanálise com o RNA publicado

19 setembro 2026. **As nove amostras reprocessadas correspondem melhor às colunas publicadas com o mesmo rótulo em todas as seis avaliações.** Não foi encontrada evidência de troca de colunas entre a nossa reanálise e a matriz publicada. Isto não verifica identidade física dos cães nem exclui erros de origem partilhados pelos dois conjuntos.

## Dados e método

Comparada a matriz não orientada corrigida v1.22.7 com GSE247355_Car_T-Dual_vs_Single_Normalized_counts_CPM.txt.gz. Dos 14.385 símbolos publicados, 11.976 têm correspondência exata e única no campo gene_name do índice ROS. Os 2.409 restantes foram excluídos desta comparação, com razões guardadas; isto não significa genes ausentes no animal. Não foram forçadas equivalências por aliases nem afirmada identidade das definições exónicas entre assemblies.

CPM reprocessado usa como denominador a soma de pares atribuídos a genes. Compararam-se todos os 81 pares de colunas: Spearman da abundância e Pearson de log2(CPM+1) após subtrair, em cada gene e separadamente em cada matriz, a média das nove amostras. A segunda análise reduz o efeito das diferenças de abundância entre genes e avalia melhor a correspondência dos perfis entre amostras. O +1 serve exclusivamente para esta transformação descritiva, sem teste diferencial.

Repetiu-se em genes com CPM≥1 em pelo menos três amostras de cada matriz (9.690 genes) e CPM≥5 (8.357 genes). São análises de sensibilidade da concordância, não novos filtros de candidatos. A atribuição global um-para-um que maximiza a soma de correlações também preserva todos os rótulos nas seis avaliações.

| Genes comparados | Métrica | Correlação entre rótulos iguais | Melhor correspondência individual |
|---|---|---:|---:|
| 11976 (all_exact) | spearman_abundance | 0.9336–0.9437 | 9/9 |
| 11976 (all_exact) | pearson_gene_centered_log2_CPM_plus1 | 0.9005–0.9548 | 9/9 |
| 9690 (CPM1_both_at_least3_samples) | spearman_abundance | 0.9097–0.9268 | 9/9 |
| 9690 (CPM1_both_at_least3_samples) | pearson_gene_centered_log2_CPM_plus1 | 0.9226–0.9698 | 9/9 |
| 8357 (CPM5_both_at_least3_samples) | spearman_abundance | 0.8875–0.9103 | 9/9 |
| 8357 (CPM5_both_at_least3_samples) | pearson_gene_centered_log2_CPM_plus1 | 0.9370–0.9801 | 9/9 |


## O que fica esclarecido

A ligação operacional entre runs, nomes de amostras da reanálise e colunas publicadas tem suporte empírico. Os resultados são compatíveis com diferenças esperadas de anotação e método de contagem; a correlação elevada não prova equivalência quantitativa gene a gene. Esta comparação usa os mesmos dados biológicos, portanto não é replicação independente.

## O que permanece aberto

Os quatro conflitos source/tissue do GEO continuam preservados em metadata_conflicts_preserved.tsv. Não se alteraram rótulos nem se declarou resolvida a identidade física. A confirmação por genótipos ou esclarecimento do depositante seria evidência adicional diferente desta concordância; não foi enviado contacto.

O protocolo declarado SMART-Seq v4 3’ DE, a história de processamento dos FASTQs e a orientação continuam a exigir reconciliação. A consistência encontrada não justifica inventar trimming de índices/UMIs nem transformar o equilíbrio forward/reverse em confirmação definitiva do kit.

Antes de inferência: documentar normalização, correspondência génica e desenho pareado de três cães. As nove bibliotecas não são nove animais independentes. Não foram calculados p-values nem declarada ausência de efeito pela simples proximidade de expressão.

Para safe harbor: a evidência RNA fornece contexto em CAR-T; continua a faltar demonstração de acessibilidade e desempenho/segurança da integração nos intervalos candidatos. Nenhum locus foi promovido ou excluído apenas por esta análise.

## Artefactos

concordance_summary.json; cross_sample_concordance.tsv (486 correlações); exact_symbol_mapping.tsv; unmapped_published_symbols.tsv; metadata_conflicts_preserved.tsv. Scripts e hashes preservados nesta pasta. Última recontagem validada: v1.22.7. Não foi iniciado novo download nem alinhamento; acompanhamento da recontagem permanece pausado porque terminou.
