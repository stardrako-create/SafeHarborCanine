# v1.21.8 — reconstrução ATAC e reavaliação verificadas

Fecho desta etapa técnica em 2026-09-14T14:41:45. Não é fecho da validação biológica do projeto.

## Resultado

A reconstrução global e o rescoring foram concluídos e revistos. **As mesmas três regiões passam o scorer legado; nenhuma é um safe harbor CAR-T validado.** Na avaliação integrada, as três regiões completas têm anotações regulatórias que exigem resolução local. A janela específica de LOC não tem sobreposição EpiC, mas a região que a contém tem.

| Região | Score legado atualizado | Percentil ATAC da região | Percentil da janela local anterior de 1 kb | Promotor/enhancer EpiC na região completa |
|---|---:|---:|---:|---:|
| ANO2/NTF3 | 0,5073 | 72,70 | 93,98 | 9.616 bp |
| LOC119876429/LOC119872513 | 0,3984 | 57,43 | 96,03 | 5.349 bp |
| NPNT/TBCK | 0,5015 | 57,63 | 96,57 | 14.800 bp |

Percentis regionais e locais respondem a perguntas diferentes, com os respetivos fundos de comprimento correspondente. Janelas escolhidas por máximos anteriores continuam exploratórias. Estes números não medem percentagem de cromatina aberta e não validam células T/CAR-T.

## Verificação da reconstrução

76 cães, 376 sequências e 2.396.858.295 bp, em bins de 25 bp. Foram verificadas as invariantes de **95.874.515 bins** nos caches e as coberturas dos BigWigs finais. Fundo de ganho reproduzido: 0,08076513558626175. Média observada em 2.328.830.335 bp e variabilidade com ≥2 cães em 2.312.901.198 bp. Ausência de cães com evidência mantém-se gap, não zero medido.

Quatro blocos foram recuperados de raw BigWigs. A conversão raw→CPM foi conferida em múltiplas sub-regiões legíveis de cada bloco; as sub-regiões ainda ilegíveis no CPM foram registadas e não usadas como falsa confirmação. Em 20 intervalos independentes (815 bins), o cálculo direto por base de média, IQR e número de cães concorda com os novos tracks, incluindo as três janelas e regiões de recuperação. Isto verifica os cálculos e escalas, não certifica os alinhamentos/experimentos originais.

**18 testes de regressão passaram**, incluindo um teste que escreve e volta a ler BigWigs com zero, um e dois cães: média ausente com zero, variabilidade ausente com menos de dois, e valores corretos quando há evidência. Corrigidos no código principal: abortar perante bloco irrecuperável em vez de usar zero; explicitar componentes de score ausentes; excluir bases de padding da média do último bin; preservar gaps na escrita de médias. O construtor rápido usado nesta execução e o código principal estão preservados. Os seus caminhos de leitura diferem; a equivalência numérica foi verificada nos intervalos de teste, não por reexecutar o construtor lento integralmente.

## Comparação com a versão anterior

461 regiões: 457 excluídas, uma com evidência insuficiente e três passes legados. **Zero mudanças no estado final ou no hard veto agregado.** Houve mudança em 205 scores finais, 246 componentes de estabilidade ATAC e 13 decisões do filtro isolado de acessibilidade; essas 13 não alteraram o estado agregado porque outros critérios permanecem ativos. As diferenças completas estão em scoring_changes.tsv. A nova regra de componentes ausentes não alterou o estado final destes 461 registos, mas fecha um comportamento incorreto para entradas incompletas.

Foi criado candidates_integrated_review_v1218.tsv: conserva o resultado legado e acrescenta a revisão regulatória regional e ausência de validação funcional. Não transforma um PASS de score em aprovação de safe harbor. Os conflitos regionais não eliminam automaticamente todos os subintervalos possíveis.

## Alternativas locais

As 2.091 janelas foram atualizadas com os novos tracks, cobertura e contagem de cães. As três janelas selecionadas anteriormente conservam cobertura ATAC observada de 100%. Cães por base/bin, mínimo/mediana/máximo: ANO2 26/33/42; LOC 16/34/57; NPNT 12/28/62. Isto é presença acima do gate, não independência entre bins nem medida validada de atividade em T.

Para orientar trabalho futuro, há 306 janelas de 1 kb sem sobreposição regulatória registada no filtro EpiC + BED externo e com covariáveis disponíveis para comparação: 116 ANO2, 115 LOC, 75 NPNT. Um conjunto de 64 janelas apresenta compromissos não dominados dentro de cada região entre maior ATAC, maior cobertura RRBS, menos repetições e menor metilação observada. Nenhum peso/limiar novo foi escolhido para forçar uma passagem. Não são 64 loci independentes: muitas janelas sobrepõem-se. A maior parte não tem as verificações de variantes, remapeamento e dados T feitas nas três janelas anteriores. Não são nomeações finais para inserção. Ver local_tradeoff_frontier_exploratory.tsv.

## O que permanece aberto

1. Identidade celular/QC completos no multiome, replicação em T caninas saudáveis/ativadas e contexto CAR-T. Os resultados anteriores de 0/1/0 fragmentos no grupo restrito de 58 núcleos permanecem insuficientes; a reconstrução do sangue não altera esses dados independentes.
2. Calibração adequada da conservação e explicação completa das pequenas diferenças de reprodução. O teste de escala em segunda região é sensibilidade, não validação de um modelo neutro.
3. Confirmação da grande SV NPNT e sequência individual; os catálogos populacionais não substituem genotipagem do animal relevante.
4. Resolução dos conflitos regulatórios e avaliação completa das alternativas antes de escolher qualquer ponto de inserção.
5. Demonstração experimental de expressão estável e neutralidade funcional após integração. Não está resolvida por este pipeline.

Por isso, a etapa de reconstrução/validação/rescoring está concluída; o objetivo de demonstrar um safe harbor canino para CAR-T não está. A conclusão continua a ser shortlist exploratória com riscos explícitos, sem vencedor validado.

## Ficheiros auditáveis

ATAC/build_manifest.json; rebuild_targeted_validation.json; rebuild_global_validation.json; test_summary.json; scoring_command.json; scoring_execution.json; scoring.log; scoring_changes.tsv; scoring_comparison_summary.json; candidates_scored_v1218.tsv; candidates_integrated_review_v1218.tsv; local_windows_evidence_status.tsv; local_refresh_summary.json; local_tradeoff_frontier_exploratory.tsv. O relatório da reunião original fica preservado como retrato anterior; o relatório atualizado contém o estado agora verificado.
