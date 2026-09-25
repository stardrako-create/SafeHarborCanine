# v1.21.4 — avaliação local, reconciliação CpG e auditoria de conservação

Data: 2026-09-13. Objetivo: candidatos para CAR-T canina. Esta versão acrescenta uma avaliação local à shortlist v1.21.3 e corrige CpG no código principal. Não é novo ranking global, reconstrução ATAC ou validação experimental. A shortlist global continua a ser a v1.21.3: três passam verificações registadas e um tem conservação em falta.

## 1. Alterações concretas e verificadas

- Corrigido `scripts/build_cpg_islands.py`: chamadas GGF e TJ mantêm coordenadas e estatísticas independentes, incluindo quando se sobrepõem. Não há etiquetas duplas derivadas de sobreposição.
- Promovido para `05_SHIP/cpg_islands_ROS_Cfam_1.0.bed` o callset integral já calculado e validado na v1.21.1: 153 908 chamadas. Código e BED anteriores guardados em `prior_main/`. As versões numeradas históricas não foram modificadas.
- Revalidado todo o BED contra comprimento, GC, Obs/Exp, tipos e duplicados; três novos testes de regressão passaram. O código principal anterior emitia 5 635 etiquetas duplas cujas coordenadas falhavam TJ. A reconciliação corrige essa classe de erro.
- Avaliadas 2 091 janelas locais dentro dos três candidatos, com comprimentos predefinidos de 500, 1 000 e 2 000 bp e passo de 250 bp (inclui a última janela completa). Todos os resultados foram conservados.
- Recalculados os máximos das médias móveis phyloP de 50 bp a partir dos WIG originais: os três valores reproduzem-se. Auditada também a chain e mapeadas as bases pontuadas para ROS_Cfam_1.0.

## 2. Resultado local: acessibilidade é heterogénea dentro das regiões

O background usa 3 000 tentativas por comprimento, seed 42, cromossomas NC_ ponderados pelo número de inícios possíveis. Estatística de acessibilidade: média exata por base no agregado ATAC/full76 existente. Os percentis são relativos a janelas aleatórias do mesmo comprimento, não a promotores, a T cells purificadas ou a um background ajustado por GC/mappability/cobertura. Não são comparáveis diretamente aos percentis das janelas regionais de 50–75 kb.

Tabela: janela de 1 kb com maior ATAC em cada região. Coordenadas BED 0-based half-open, assembly ROS_Cfam_1.0.

| Candidato | Janela de 1 kb | Percentil local | Repetições | RRBS observado | Metilação nas bases observadas | Máximo de cães com pico ATAC |
|---|---|---:|---:|---:|---:|---:|
| ANO2/NTF3 | NC_051831.1:39728324–39729324 | 94.10 | 0.0% | 27.4% | 64.79% | 0 |
| LOC119876429/LOC119872513 | NC_051811.1:17235706–17236706 | 96.07 | 9.4% | 20.0% | 59.53% | 3 |
| NPNT/TBCK | NC_051836.1:27057411–27058411 | 96.63 | 13.7% | 8.9% | 100.00% | 9 |

As janelas são exemplos exploratórios selecionados depois de procurar máximos, não sítios de inserção aprovados. A pesquisa de muitos subintervalos aumenta a probabilidade de encontrar percentis altos; estes percentis não são p-valores corrigidos para a seleção. A seleção deve ser confirmada por dados independentes ou reamostragem por cão antes de se tornar critério de avanço.

Os três exemplos de 1 kb não sobrepõem CpG, BED regulatório externo, fronteiras TAD ou picos consenso nos inputs usados. Isto não demonstra ausência de função regulatória: as anotações continuam incompletas e a ausência de pico consenso coexiste com picos individuais pouco frequentes.

ANO2: a janela de 1 kb mais acessível tem 0 bp de RepeatMasker, p94,1 e RRBS em 27,4% das bases; não há pico individual nessa janela na track de frequência. O melhor máximo de 500 bp de ANO2 sobrepõe 51,6% de repetição, mostrando por que não se deve escolher só pelo ATAC.

NPNT: o melhor 1 kb atinge p96,63, mas RRBS cobre apenas 8,9%, com metilação de 100% nas bases observadas e máximo de três cães no indicador de cobertura. O melhor 500 bp não tem RRBS observado. Isto não descreve a metilação de toda a janela e não permite concluir silêncio da cassete.

LOC119876429: o melhor 1 kb atinge p96,07, com 9,4% de repetição e RRBS em 20% das bases. A média observada é 59,53%. O máximo local de cobertura RRBS chega a 60 cães; este máximo não significa 60 cães medidos em todas as bases.

Os exemplos mantêm a escala da Mother Track própria: não são CPM convencional nem demonstração de acessibilidade de CAR-T. Máximo de frequência de pico também não equivale a número de cães com sinal ATAC medido. Gaps da track esparsa de frequência são interpretados como zero, conforme a semântica usada pelo scorer.

## 3. Conservação: integridade numérica confirmada, modelo ainda piloto

| Candidato | Máximo média móvel 50 bp reproduzido | Fração da região ROS com phyloP mapeável |
|---|---:|---:|
| ANO2/NTF3 | 5,10122 | 99,84% |
| NPNT/TBCK | 4,22592 | 99,93% |
| LOC119876429/LOC119872513 | 4,01352 | 99,89% |

Cada um tem uma única chain com bases alinhadas na região consultada. Os valores WIG foram associados às coordenadas CanFam3 usando a chain; bases sem correspondência mantêm-se desconhecidas. Os máximos por base nas tabelas locais são uma estatística diferente do máximo de média de 50 bp e não devem ser comparados diretamente com o veto 6,5.

Estas verificações sustentam integridade e cobertura do cálculo existente. Não validam o modelo evolutivo: continua `neutral_model_1region.mod`, ajustado a uma região de 100 kb. Não foi treinado um novo modelo neutro nesta execução. A escolha de regiões neutras independentes, sensibilidade do modelo, suporte de espécies por base e eventual concordância com outra referência conservacional continuam necessários antes de usar o valor como evidência definitiva.

### Correção importante sobre o quarto candidato, LOC111090199/LOC100682550

A leitura direta da chain encontra 16 chains com algum alinhamento, mas a melhor é a chain 21 para chr16: **62 605 das 71 869 bases da região ROS estão alinhadas (87,11%)**, em 54 blocos. O intervalo-alvo abrangido mede 63 693 bp, dos quais 98,29% correspondem a bases alinhadas. É o mesmo intervalo chr16:21250606–21314299 registado no liftOver permissivo.

Portanto, não é correto concluir, a partir dos resultados anteriores, que só existem fragmentos menores de 1 kb em cromossomas dispersos ou que está demonstrada uma discordância estrutural biológica. Existe uma correspondência dominante substancial, coexistindo com alinhamentos menores. Os ~12,89% não alinhados na fonte e os blocos/gaps impedem tratá-la como equivalência integral sem avaliação adicional.

O candidato permanece sem conservação utilizável no scoring: não calculámos phyloP para ele nem promovemos uma aproximação a evidência completa. Próximo passo justificado: verificar sequência, contexto e reciprocidade da correspondência dominante, e depois avaliar conservação nas bases efetivamente alinhadas com cobertura explícita. A causa dos gaps não foi determinada.

## 4. Limites que não desapareceram

- PBMC/full76 não demonstra o estado cromatínico de células T caninas ativadas/CAR-T. Não foi feita análise por cão nesta versão.
- O veto ATAC na janela regional permanece no scorer. A avaliação local não mudou filtros para recuperar nomes; documenta o sinal dentro da shortlist atual.
- Os picos consenso usados são os do input legado de 71 cães, enquanto a acessibilidade/frequência usam full76. Esta proveniência está explícita no manifest anterior e deve ser reconciliada antes de congelar o protocolo.
- Repetições locais foram obtidas dos outputs RepeatMasker reais, com intervalos 1-based inclusivos convertidos para BED e união das sobreposições. Isso não substitui uma avaliação local de mappability, variantes, estrutura e genes/transcritos.
- A janela regional passar as distâncias não equivale a uma revisão independente de todas as anotações locais. Não foram desenhados gRNAs nem escolhido um ponto de integração.
- O modelo phyloP continua piloto; os dados ATAC continuam anteriores à reconstrução completa de blocos recuperados.

## 5. Ficheiros, validação e próximos passos

`all_local_windows.tsv` preserva as 2 091 janelas; `top_local_windows_exploratory.tsv` mostra três máximos não sobrepostos por candidato/escala, sem declarar segurança. `local_backgrounds.tsv` permite interpretar os percentis. `conservation_audit.tsv`, `chain_alignment_audit.tsv` e `chain_blocks.json` documentam a conservação. `cpg_reconciliation.json` regista a promoção do BED e o backup. `validation.json` documenta testes e verificações.

Prioridades seguintes: confirmação de acessibilidade por cão nas janelas locais e coerência 71/76; avaliação local de mappability/variantes e anotação regulatória; validação independente do modelo de conservação; análise da correspondência dominante do quarto locus; confirmação em contexto celular CAR-T. A escolha entre ANO2 e NPNT não ficou resolvida pelo score global, nem pelo máximo local.
