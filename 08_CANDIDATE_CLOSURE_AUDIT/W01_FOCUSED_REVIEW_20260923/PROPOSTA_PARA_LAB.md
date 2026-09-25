# Proposta para discussão: w01 e dois comparadores

23 de setembro de 2026 · revisão computacional dirigida

## Decisão proposta

Levar **w01 como primeiro locus a testar**. Levar **bg1k_2848 como primeiro comparador genómico a discutir**, com **bg1k_0442 como segundo comparador ou alternativa**, conforme a viabilidade do desenho nos cães escolhidos. Manter w11 condicionado à resolução da dúvida estrutural. Não é necessário resolver w11 para preparar o trabalho de w01.

Esta ordem é uma proposta de planeamento, sem novo limiar de aprovação. O 2848 combina ausência de repetidos detetados, cobertura regulatória/conservação completa nas análises feitas e menor diferença de GC face a w01. O 0442 oferece uma segunda região, mais distante do início do transcrito codificante mais próximo, mas contém 128 bp de repetidos e maior diferença de GC. A escolha pode mudar quando houver nuclease, cassete e sequência individual.

## O que sabemos dos três loci

Todas as janelas têm 1.000 bp. Coordenadas ROS_Cfam_1.0 em BED: início zero, fim exclusivo.

| Medição | w01 | bg1k_2848 | bg1k_0442 |
|---|---:|---:|---:|
| Cromossoma | NC_051811.1 | NC_051840.1 | NC_051837.1 |
| Intervalo BED | 17255706–17256706 | 28808728–28809728 | 10430486–10431486 |
| Percentil ATAC no sangue | 63,75 | 15,21 | 7,03 |
| GC | 41,1% | 35,9% | 31,7% |
| Repetidos detetados | 0 bp | 0 bp | 128 bp |
| Distância ao início anotado de RNA mais próximo | 32.785 bp | 82.639 bp | 69.671 bp |
| Distância ao início do transcrito codificante mais próximo | 110.509 bp | 140.510 bp | 710.285 bp |
| Alelos normalizados / PASS / PASS com AF≥1% | 27 / 25 / 10 | 24 / 20 / 9 | 35 / 30 / 13 |
| Segmentos de 20 bp únicos, sem sobrepor REF de variantes PASS | 567 | 579 | 534 |
| Sobreposições no rastreio por intervalos do catálogo SV | 0 | 0 | 0 |

As contagens de variantes são entradas do catálogo, não genótipos dos cães nem necessariamente variantes independentes entre fontes. Os segmentos de 20 bp são possibilidades de sequência para avaliação posterior; ainda faltam PAM, mismatches, bulges e especificidade dependente da nuclease.

Os três loci têm correspondência exata integral ROS–UU na análise usada para variantes. Não há sobreposição de transcritos nesta revisão. O RNA mais próximo é lncRNA nos três casos. As distâncias são à extremidade 5′ anotada: inícios alternativos em células T podem estar incompletamente anotados.

O contexto TAD não teve genes de risco catalogados nos proxies testados às várias escalas. Estes proxies não medem contactos cromatínicos nas células T dos cães. Nas regiões EpiC consultadas não foram encontrados estados 1–7, usados no projeto como promotor/enhancer. O 2848 apresenta estado 9; não deve ser descrito como exclusivamente quiescente. O 0442 apresenta estado 13 nas fontes revistas. Ambos têm 1.000 bp de scores de conservação com suporte de outra espécie; cobertura completa não prova neutralidade.

## A dúvida principal de w01

A evidência de acessibilidade em células T continua escassa: dois inícios de fragmentos em dois núcleos T, num único dador. O percentil de sangue não substitui essa medição. RRBS observou apenas 20% da janela, com metilação média de 60,71% nas bases observadas. O lncRNA LOC119872513 começa, pela anotação, a 32.785 bp; o seu corpo fica a 17.805 bp. Importa avaliar o contexto regulatório local, para além da distância ao gene codificante.

É isto que mais pode alterar a escolha de w01: atividade no material celular relevante, diferenças individuais que afetem o desenho e efeitos sobre os transcritos vizinhos. As pendências dos restantes controlos não são automaticamente problemas deste locus.

## Para que servem os comparadores

2848 e 0442 vieram do conjunto aleatório condicionado pelos filtros iniciais. Têm ATAC de sangue baixo, mas não foram estabelecidos como controlos negativos em células T. A comparação com w01 pode explorar diferenças de comportamento entre loci; três loci com GC e contexto distintos não isolam um efeito causal da acessibilidade.

Estes comparadores genómicos complementam os controlos do procedimento. Manter células do mesmo cão sem intervenção e um controlo da entrega adequado ao sistema escolhido. Se o laboratório já tiver um comparador de integração caracterizado, discutir a sua inclusão. O cão é a unidade biológica; clones e medições técnicas devem ficar identificados dentro de cada cão.

## Decisões a fechar na reunião

| Decisão | Informação necessária | Consequência |
|---|---|---|
| Sistema de edição | Nuclease/variante e estratégia de integração/entrega | Permite verificar alvos e especificidade |
| Cassete | Sequência completa, versão e objetivo de expressão | Define o desenho e o que medir |
| Material biológico | Cães/amostras e dados de sequência disponíveis | Permite confirmar cada locus e variantes relevantes |
| Comparadores | Um ou dois loci adicionais e controlos já usados pelo laboratório | Define a comparação que a experiência pode sustentar |
| Critérios de sucesso | Integridade/cópias, expressão persistente e heterogénea, genes vizinhos, viabilidade e função | Fixar critérios antes de interpretar os resultados |

Uma alteração estrutural individual, um alvo inviável para a nuclease escolhida ou um efeito regulatório relevante obriga a reconsiderar o locus correspondente. Não foram escolhidos guias, primers ou condições de edição nesta revisão.

## Evidência e limites da revisão

Os dez hashes das fontes da revisão integrada foram reconferidos. O filtro de segmentos únicos/PASS foi aplicado aos dois controlos. O rastreio SV usou sobreposição de intervalos no VCF local Dog10K; não incluiu análise de mates BND ou intervalos de confiança dos breakpoints, e ausência de registos não exclui SV individual. Foi relida a anotação completa para comparar os inícios de todos os tipos de RNA nos três loci.

Ficheiros nesta pasta: focused_checks.json, all_RNA_context.json e source_manifest.json. As análises anteriores mantêm-se intactas. A validação da integração, expressão e função requer trabalho experimental.

## Complemento SV — 24/09/2026
O varrimento integral de 141.810 registos incluiu as posições remotas das 150 ALT BND. Não encontrou sobreposições ou ligações BND às janelas w01, 2848 e 0442. Recuperou o registo de deleção já conhecido em w11, cuja interpretação permanece aberta. A fonte não fornece CIPOS/CIEND. Recibos: SV_extended_20260924/status.json e validation.json. Este complemento amplia a consulta por intervalos descrita acima; não exclui variantes estruturais individuais.
