# v1.23.1 — teste dos mates concluído e revisto

19 setembro 2026. **27/27 alinhamentos concluídos e verificados independentemente. R1 e R2 fornecem alinhamento genómico e atribuição génica substanciais nas nove amostras.** O teste não sustenta descartar R2 como se fosse composto apenas por índices/poly(T) nos ficheiros tratados analisados.

## Desenho e validação

Por amostra, os mesmos primeiros 100.000 pares tratados foram alinhados em três modos: paired, apenas R1 e apenas R2. Nove amostras; 900.000 pares distintos entre as amostras, reutilizados nos três modos. Não são 2,7 milhões de pares independentes. Índice corrigido v1.22.7, STAR 2.7.11b. Não foram alterados parâmetros do pipeline completo nem contagens v1.22.7.

O revisor verificou os 27 hashes de contagens, hashes dos subsets contra recibos, finalização STAR, 42.309 IDs esperados, soma das categorias GeneCounts=100.000 unidades em cada orientação, concordância das tabelas com os logs e recalculou as 27 correlações. Todas as verificações passaram. A unidade single-end é um read; no modo paired é um par. As taxas abaixo usam esses denominadores; não somar contagens R1+R2 como moléculas independentes.

## Resultados nas nove amostras

| Modo | Mapeamento único (%) | Atribuição génica (%) |
|---|---:|---:|
| paired | 88.99–93.19 | 74.810–79.506 |
| R1 | 87.63–91.62 | 73.431–78.120 |
| R2 | 87.75–91.88 | 73.535–78.106 |


Correlações de Spearman de contagens não orientadas, restringindo cada comparação a genes não zero em pelo menos um dos modos:

- R1_paired: 0.9754–0.9828.
- R2_paired: 0.9751–0.9812.
- R1_R2: 0.9657–0.9749.

Os pares têm maior taxa de mapeamento único e atribuição génica que cada mate isolado neste teste. A proximidade de R1 e R2 reforça que ambos contêm sequência utilizável. Não há fundamento neste resultado para remover R2 nem cortar seis bases automaticamente. A concordância é parcialmente esperada porque os mates vêm dos mesmos fragmentos e não constitui replicação biológica independente.

## O que fica por resolver

O teste fornece suporte empírico ao uso operacional dos pares disponibilizados, mas não identifica o kit original, não reconstrói desmultiplexagem/processamento anterior e não resolve os quatro conflitos de metadados. A discrepância com a estrutura descrita no manual continua documentada na v1.23.0. O presente teste usou reads tratados; a auditoria anterior avaliou a estrutura dos brutos disponíveis.

Amostragem inicial determinística, não aleatória, e profundidade de apenas 100.000 pares por amostra. Não interpretar zeros de genes raros, diferenças pequenas entre modos ou resultados dos loci candidatos como testes de estabilidade, equivalência ou segurança. Não foram gerados p-values. A contagem completa v1.22.7 e o contexto dos 40 genes v1.22.9 permanecem as referências para descrição de expressão.

Não foi demonstrada acessibilidade em CAR-T dos intervalos candidatos, persistência do transgene, ausência de perturbação regulatória nem segurança da integração. Nenhum candidato foi promovido a safe harbor.

## Espaço e fecho

Pasta desta análise: 0.713 GB, abaixo da estimativa de 1 GB. Espaço livre observado no fecho: D 129.0 GB e C 61.1 GB. A queda total de espaço em D durante o intervalo excede o tamanho desta pasta; a origem do restante consumo não foi determinada e não foi atribuída a este teste. Não foram apagados ficheiros de outros processos nem inputs não validados. O worker terminou.

O acompanhamento será pausado após este fecho. Artefactos: independent_review.json, execution_complete.json, mate_alignment_QC.tsv, mate_count_concordance.tsv e recibos por amostra/modo. Scripts mate_sensitivity_v1231.py, review_mates_v1231.py e finish_mates_v1231.py.
