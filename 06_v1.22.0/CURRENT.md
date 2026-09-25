# Safe Harbor CAR-T — v1.22.0: ATAC nas alternativas locais

14 de setembro de 2026. **A análise foi alargada às 64 janelas alternativas. Há sinais T esparsos fora das três janelas anteriores, mas continuam a faltar suporte replicado e validação funcional. Não há vencedor demonstrado.** O scoring global continua a ser o da v1.21.8; esta versão acrescenta evidência local.

## O que foi feito

Mapeámos as 64 janelas de 1 kb por duas vias: ROS→CanFam3→CanFam6 e inversão da cadeia CanFam6→ROS. Exigimos concordância por base e unicidade do destino. Em 37 janelas, o mapeamento é integral e contínuo. Nas outras 27, há pequenas lacunas, diferenças entre vias ou descontinuidades: contamos separadamente apenas inícios em bases de destino com concordância, guardando explicitamente a cobertura. Não imputámos zero às bases não avaliáveis.

Foram consultados três intervalos regionais de fragmentos do [multiome público GSE244116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244116), reutilizando o cache. As consultas foram bracketed por registos ordenados antes/depois dos intervalos. A descompressão independente com HTSlib reproduziu os 2.027, 4.613 e 5.167 registos guardados para os intervalos de LOC, ANO2 e NPNT. A completude das contagens de início depende da ordenação do ficheiro; não é uma contagem exaustiva de fragmentos longos sobrepostos cujo início esteja fora da consulta.

Usámos os mesmos grupos definidos por RNA na v1.21.9: 155 núcleos na interseção das quatro análises, 157 na união, 5.849 núcleos QC no total. Em 61 das 64 janelas, as contagens são iguais entre interseção e união. Em três janelas sobrepostas ANO2 (inícios ROS 39729824, 39730074 e 39730324), a união acrescenta um início por janela onde a interseção tem zero. Assim, são 13 janelas positivas na interseção e 16 na união; a tabela abaixo refere-se à interseção de 155. São núcleos de um único tumor, não réplicas biológicas ou T caninas saudáveis validadas.

## Resultado

| Região | Alternativas | Mapeamento integral/contínuo | Menor cobertura de concordância nas alternativas | Janelas com ≥1 início no grupo T, contando as bases concordantes | Fragmentos distintos / núcleos distintos em todas as alternativas da região |
|---|---:|---:|---:|---:|---:|
| ANO2/NTF3 | 22 | 13 | 99,7% | 1 | 1 / 1 |
| LOC119876429/LOC119872513 | 12 | 5 | 92,3% | 3 | 3 / 3 |
| NPNT/TBCK | 30 | 19 | 93,6% | 9 | 9 / 9 |

As 13 janelas com algum sinal não são 13 descobertas independentes: há sobreposição e reutilização dos mesmos fragmentos entre janelas. Os totais da última coluna removem essa duplicação dentro de cada região. As bases sem correspondência continuam sem avaliação; contagens nas bases concordantes não equivalem necessariamente a uma janela integral de 1 kb.

Entre as 37 janelas integralmente mapeadas, têm algum sinal zero de 13 em ANO2, duas de cinco em LOC e oito de 19 em NPNT. A inclusão cuidadosa das bases concordantes das restantes 27 acrescenta uma janela positiva em cada região.

## Janelas úteis para decidir a próxima verificação

Estas são exemplos de compromissos, não locais aprovados para edição. Coordenadas ROS, 0-based half-open. Todas pertencem à frente exploratória anteriormente definida sem usar estas contagens T; priorizá-las agora com estes mesmos dados continua a exigir validação independente.

| Janela | Inícios no grupo T de 155 | Inícios em todos os 5.849 QC | Percentil ATAC no sangue | Fração RRBS observada | Máximo da média phyloP publicada em 50 bp |
|---|---:|---:|---:|---:|---:|
| LOC: NC_051811.1:17255706–17256706 | 2 | 18 | 63,75 | 20,0% | 0,342 |
| NPNT: NC_051836.1:27041161–27042161 | 2 | 24 | 89,50 | 10,0% | 0,314 |
| NPNT: NC_051836.1:27050411–27051411 | 4 | 490 | 66,61 | 3,9% | 1,165 |

Estas três janelas têm mapeamento integral concordante. A janela NPNT com quatro inícios T tem também 490 inícios no conjunto QC: quatro não demonstra enriquecimento específico em T. Não foi feita normalização pela profundidade genómica de fragmentos nem inferência estatística de enriquecimento. A cobertura RRBS de 3,9% é particularmente limitada; o resto não pode ser chamado desmetilado.

A primeira janela LOC oferece um ponto concreto para aprofundar a revisão local pela menor conservação e ausência de sobreposição regulatória no filtro que definiu a frente. Duas moléculas em dois núcleos continuam a ser suporte muito escasso. A janela NPNT com percentil sanguíneo 89,50 permite estudar outro compromisso, também com só dois núcleos e RRBS incompleto. Nenhuma está sujeita ainda a todas as verificações individuais de variantes/mapeabilidade efetuadas nas três janelas originais. A chamada SV regional NPNT mantém a interpretação revista na v1.21.9: catalogada e questionada por heterozigotia distribuída, ainda por resolver.

## O que isto muda

Os zeros das três janelas anteriores não se generalizam a todas as alternativas. Há sinais focais noutros subintervalos, mas a dimensão da evidência continua insuficiente para afirmar acessibilidade robusta em T ou segurança para CAR-T. Procurar apenas o máximo de ATAC no sangue não seleciona necessariamente a janela com melhor suporte na população T exploratória.

A próxima etapa computacional útil é aplicar às alternativas priorizadas a revisão de variantes locais, mapeabilidade e contexto regulatório com a mesma exigência das janelas originais, e completar o QC/identidade do multiome. Para fechar a conclusão biológica, continuam indispensáveis dados de T caninas relevantes com replicação e validação da integração/neutralidade funcional. Não se deve baixar filtros para promover estes sinais escassos a candidatos validados.

## Verificação e ficheiros

Foram verificadas 64 linhas, limites de cobertura, inclusão das contagens da interseção nas da união e identificação explícita das três diferenças e equivalência das contagens por bases concordantes com as contagens integrais nas 37 janelas aplicáveis. A descompressão HTSlib coincide com os registos consultados. Os critérios anteriores de seleção e scores foram preservados.

`frontier_T_evidence_v1220.tsv` contém todas as 64 janelas, cobertura e contagens integrais/parciais separadas; `positive_windows_exploratory.tsv`, as 13 com algum sinal; `frontier_shared_mapping_blocks.json`, os blocos auditáveis; `frontier_mapping_summary.json`, o critério de mapeamento integral; `frontier_T_summary.json`, a análise que acrescenta as bases concordantes das restantes janelas; `regional_fragment_queries.json`, a proveniência das consultas. Os intervalos com concordância parcial não receberam aprovação de mapeamento integral.
