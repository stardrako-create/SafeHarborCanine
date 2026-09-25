# Recontagem corrigida v1.22.7 — conclusão e revisão, 19 setembro 2026

**Execução computacional concluída e consistência dos outputs verificada nas nove amostras. A validação biológica de safe harbor continua por fazer.** Este estado substitui os checkpoints anteriores de execução.

## Correção e validação

Corrigido exclusivamente o campo source do GTF que continha espaços (`Curated Genomic`), incompatível com a tokenização do STAR. Coordenadas e identificadores mantidos. O novo índice contém os 42.309 genes esperados, incluindo os 60 antes omitidos (61 exões). Referência ROS_Cfam_1.0 e STAR 2.7.11b; esta não é reprodução exata do pipeline original CanFam3.1/HTSeq.

Recontadas as nove amostras, sequencialmente, reutilizando os pares tratados já disponíveis. O worker verificou hashes de referência, índice e inputs tratados, finalização STAR e coerência entre reads pós-fastp, pares de entrada e soma das categorias GeneCounts. Uma revisão separada verificou novamente hashes dos outputs e matrizes, 42.309 IDs por matriz e igualdade célula a célula com cada ficheiro individual nas três orientações. Os certificados MD5 dos downloads tinham sido auditados na v1.22.6; não se alegou novo MD5 integral dos brutos nesta revisão.

## Impacto medido

Os 60 genes recuperados recebem **140,423 pares atribuídos no total das nove amostras** (contagens não orientadas). Isto não representa novos reads sequenciados: são atribuições possibilitadas pela anotação corrigida. Entre os genes já presentes, mudaram 14 genes distintos, em 38 combinações gene/amostra; soma de diferenças absolutas = 469 contagens. As categorias não atribuídas podem também variar.

| Amostra | Pares de entrada | Mapeamento único | Pares nos 60 genes recuperados | Genes anteriores alterados |
|---|---:|---:|---:|---:|
| B_B7H3_CAR_T | 71,312,581 | 92.35% | 17,021 | 2 |
| B_BC_CAR_T | 73,105,889 | 92.08% | 21,723 | 5 |
| B_T_cell | 63,103,725 | 91.46% | 17,555 | 4 |
| E_B7H3_CAR_T | 63,968,677 | 93.26% | 18,594 | 4 |
| E_BC_CAR_T | 50,056,953 | 93.10% | 11,736 | 2 |
| E_T_cell | 47,628,160 | 93.54% | 14,694 | 5 |
| M_B7H3_CAR_T | 63,529,252 | 89.13% | 13,388 | 3 |
| M_BC_CAR_T | 58,110,849 | 90.91% | 13,595 | 8 |
| M_T_cell | 59,387,857 | 89.52% | 12,117 | 5 |


Os dez genes de contexto revistos têm contagens brutas não orientadas inalteradas em todas as amostras. A correção não autoriza promover os candidatos a safe harbors nem altera diretamente os resultados ATAC. O denominador CPM mudou, pelo que contagem bruta igual não implica CPM exatamente igual.

## Contexto génico corrigido

| Gene | Intervalo de contagens brutas nas 9 amostras | Intervalo CPM* |
|---|---:|---:|
| ANO2 | 0–4 | 0.0000–0.0702 |
| NTF3 | 0–4 | 0.0000–0.0808 |
| NPNT | 0–13 | 0.0000–0.2886 |
| TBCK | 741–1187 | 18.1774–23.6933 |
| TSEN15 | 1642–2746 | 37.1203–48.2105 |
| C7H1orf21 | 755–1823 | 16.7493–33.1226 |
| CD247 | 6800–16027 | 150.8550–322.9327 |
| CD8B | 469–17105 | 10.4046–305.6224 |
| ARPC5 | 18112–30295 | 444.2237–579.3297 |
| TET2 | 2371–3494 | 55.7942–67.0142 |


*CPM = contagem não orientada × 1 milhão / soma de contagens génicas não orientadas da mesma amostra. Não assumir igualdade com o denominador da tabela publicada. Os intervalos combinam condições e cães; não são testes de diferenças, equivalência ou estabilidade.

CD247 e CD8B podem agora ser avaliados na anotação ROS e apresentam contagens substanciais; a ausência na tabela CPM publicada não significava ausência de expressão. ANO2, NTF3 e NPNT têm contagens baixas neste conjunto. Isso não demonstra silêncio absoluto, acessibilidade da janela candidata ou segurança de integração. Genes vizinhos expressos, como TBCK e TSEN15, continuam relevantes para avaliar perturbação regulatória, mas expressão por si não demonstra interação com o intervalo.

## Limitações e próximos passos científicos

- Resolver ou manter explicitamente os quatro conflitos nos campos GEO de identidade/condição. Rótulos de biblioteca, título e desenho publicado não são uma confirmação física da identidade. Não efetuar trocas silenciosas.
- A contagem forward/reverse aproximadamente equilibrada é consistente com biblioteca não orientada; reconciliar com o kit SMART-Seq v4 3’ DE declarado e a história de processamento. Não introduzir cortes de reads só pelo nome do kit.
- Q30=100% é pouco informativo neste material com scores constantes de 30 na amostra inspecionada. Não apresentar como qualidade perfeita.
- Comparar expressão com o ficheiro publicado mediante correspondência inequívoca de genes e normalização documentada; não tratar discrepâncias ROS/CanFam3.1 ou STAR/HTSeq como erro automaticamente.
- Qualquer inferência deve respeitar três cães pareados em três condições. Não executar testes de expressão diferencial sobre CPM como se fossem contagens brutas; esta revisão é descritiva.
- RNA génico de CAR-T oferece contexto celular relevante, mas não substitui acessibilidade do locus, estabilidade de expressão do transgene, impacto em genes vizinhos e segurança genómica/funcional. Nenhum locus está validado como safe harbor.

## Espaço e proveniência

Espaço livre no fecho: D 146.0 GB; C 60.3 GB. Sem novos downloads; não foi necessário apagar inputs. Ficheiros originais e outputs anteriores preservados para reprodução. O worker terminou e não deve ser reiniciado para repetir esta execução. O acompanhamento desta recontagem será pausado após guardar este relatório.

Evidência: execution_complete.json; index_validated.json; samples/*/validated.json; independent_review.json; comparison_with_v1226.json; existing_gene_count_changes.tsv; review_sample_QC.tsv; review_focus_gene_context.tsv; gene_counts_unstranded.tsv, forward e reverse. Scripts corrected_rna_v1227.py, review_corrected_v1227.py e finish_corrected_v1227.py.
