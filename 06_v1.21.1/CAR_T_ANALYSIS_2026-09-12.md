---
title: "Safe Harbor CAR-T canina — análise integrada e critérios de fecho"
date: 2026-09-12
status: protocolo-em-consolidacao
tags: [canino, CAR-T, safe-harbor, auditoria, bioinformatica]
version_analyzed: v1.21.1
---

# Safe Harbor CAR-T canina — análise integrada e critérios de fecho

## 1. Decisão de âmbito e conclusão atual

O utilizador confirmou em 2026-09-12 que o destino é CAR-T canina e pediu um protocolo individual, independente de colaboração externa. O trabalho computacional deve identificar locais candidatos para integração de uma cassete CAR em células T caninas e documentar acessibilidade, riscos anotados e incerteza. O objetivo celular deixa de estar em aberto. Continuam por especificar composição CD4/CD8, estado de ativação, cassete/promotor e requisitos funcionais de expressão. Estes elementos não são presumidos nesta análise.

Esta nota atualiza o estado do projeto: referências antigas a V9-B2, 26 sobreviventes, confirmação independente ou shortlist robusta não descrevem o resultado atual. A v1.21.1 avaliou 461 regiões: 457 excluídas, três com evidência insuficiente e uma que passa as verificações registadas. Nenhuma é um safe harbor validado. A execução não foi uma reconstrução completa desde FASTQ.

O protocolo pode ser fechado como método reprodutível de descoberta de candidatos, com critérios e limitações explícitos. Segurança absoluta não é uma propriedade demonstrável pelo score. A validação funcional pertence a uma fase distinta e deve limitar a conclusão ao locus, cassete e contexto celular efetivamente estudados.

## 2. O que os 71 cães fornecem e o que não fornecem

Jin et al. (2024), DOI 10.1111/acel.14079, estudaram metilação e acessibilidade em PBMC de 71 cães. Portanto, a descrição precisa é PBMC de sangue periférico; não sangue total e não células T purificadas. O objetivo original era estudar envelhecimento epigenético. O artigo considera composição celular e associações com idade. Fonte: [Jin 2024](https://pubmed.ncbi.nlm.nih.gov/38263575/) e [artigo](https://onlinelibrary.wiley.com/doi/10.1111/acel.14079).

Interpretação para este projeto: a coorte é útil para descoberta e avaliação de heterogeneidade em células mononucleares sanguíneas. O sinal agregado não identifica, por si só, a população celular que o produz. Um pico forte pode ser sustentado por outra população de PBMC; um sinal fraco pode refletir diluição de sinal específico de células T. Estes são mecanismos plausíveis, não causas já demonstradas para os nossos candidatos.

Ter 71 animais não equivale a ter 71 medições informativas por locus, nem a replicação em células T ativadas ou CAR-T. Deve ser exportado o número efetivo de cães com evidência em cada medida. A mistura de células, idade, composição racial e qualidade técnica podem influenciar média e variabilidade; a importância de cada fator nos nossos loci ainda precisa de ser medida.

Há uma distinção adicional essencial: o manifest da v1.21.1 aponta para ATAC/full76, enquanto RRBS deriva da coorte original de 71 cães. Logo, não descrever a execução atual como exclusivamente baseada nos 71. A independência do protocolo não implica apagar a proveniência de dados previamente incorporados. Proponho a coorte original como análise principal predefinida e as cinco amostras adicionais como análise separada de generalização, depois de conferir metadados e compatibilidade. Isto é uma proposta; não foi executado neste turno.

O Hi-C provém de um animal, Mischka, com três bibliotecas destinadas à montagem genómica. O tipo celular da preparação Hi-C não está identificado na documentação consultada. Serve de contexto estrutural com limitações, não como mapa de TAD validado em células T caninas. Evidência humana, ou convertida entre assemblies, deve continuar marcada como transferida, sem equivalência presumida.

## 3. Estado reprodutível da v1.21.1

Raiz dos resultados: `D:/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/06_v1.21.1/`.

| Região candidata | Coordenadas ROS_Cfam_1.0, BED | Score anterior → corrigido | RRBS observado | Fração da janela observada | Estado |
|---|---|---:|---:|---:|---|
| ANO2/NTF3 | NC_051831.1:39677324–39739751 | 0,4309 → 0,5072 | 74,03% | 11,29% | Passa verificações registadas |
| NPNT/TBCK | NC_051836.1:27013911–27071204 | 0,5120 → 0,5015 | 81,49% | 8,38% | Faltam repetições e conservação |
| LOC119876429/LOC119872513 | NC_051811.1:17217706–17274511 | 0,4605 → 0,3988 | 84,71% | 6,16% | Faltam repetições e conservação |
| LOC111090199/LOC100682550 | NC_051820.1:22976673–23048542 | 0,3556 → 0,3740 | 76,74% | 7,17% | Faltam repetições e conservação |

Os scores não são probabilidades de segurança ou sucesso. A média RRBS é ponderada pelas bases dos bins observados; não representa uma média por CpG/read nem a metilação da janela inteira. O background RRBS corrigido tem mediana ~77,15%; isso explica que ANO2 melhore relativamente apesar de a metilação absoluta corrigida ser muito superior à anterior.

Correções já realizadas:

- Média RRBS condicionada à máscara de evidência, tanto em candidatos como em background; zeros medidos preservados, ausência de evidência tratada como desconhecida.
- Builder RRBS versionado passa a escrever gaps nas posições sem evidência. Tracks legadas foram interpretadas com a máscara, não reconstruídas neste run.
- Falta de anotação obrigatória deixa de dar aprovação completa; componentes de score ausentes não são compensados por redistribuição silenciosa dos pesos.
- Erros de intervalos/leitura propagam; background sem amostras suficientes termina com erro.
- Reconstrução CpG genome-wide com 112 933 chamadas GGF e 40 975 TJ, 153 908 no total. Cada chamada mantém coordenadas e estatísticas próprias. Estes são números de chamadas, não de regiões únicas não sobrepostas.
- Verificação dos limiares exportados e ausência de duplicados de coordenadas com o mesmo tipo. O algoritmo sliding/merge continua próprio do projeto; não foi demonstrada equivalência a uma implementação canónica publicada.
- Oito testes de regressão passaram, incluindo escrita real de BigWig, ausência de evidência e máscara do background.

O matching entre estatística de janela e background, o uso de leitura exata e a correção da distância aos 5′ dos genes integram a linhagem corrigida. A v1.21.1 não atribui a si a autoria de todas essas correções anteriores. Os testes efetuados não equivalem a validação científica integral.

## 4. ATAC baixo: hipóteses separadas e conflito confirmado no desenho

O código atual inclui um piso de acessibilidade p55 e, simultaneamente, um veto por qualquer sobreposição da janela com picos ATAC consenso. Também inclui uma componente que favorece baixa frequência de picos. Assim, o algoritmo favorece acessibilidade difusa suficiente, evitando picos frequentes. É um desenho próprio que pode divergir da intenção de selecionar um local acessível para CAR-T.

Uma janela de 50–75 kb é eliminada mesmo que o pico esteja longe do eventual ponto de inserção. O código confirma essa geometria; ainda falta quantificar quantas exclusões seriam diferentes numa análise local. Um pico ATAC não prova, sozinho, que o intervalo seja um enhancer ou promotor indispensável. Também não prova que seja seguro interrompê-lo. A decisão precisa de combinar evidência regulatória independente e localização concreta.

Não proponho remover indiscriminadamente este veto. Proponho medir a sua contribuição única e sobreposta, definir a unidade local de análise antes do ranking e distinguir acessibilidade de risco regulatório. Comparar uma versão predefinida local com a versão atual permite avaliar o custo do filtro sem escolher regras para recuperar favoritos.

Outras causas plausíveis do ATAC baixo: diluição pela média de janela, composição celular, normalização, universo inicial restrito e artefactos de processamento. O manifest confirma que as tracks atuais antecedem a recuperação dos blocos problemáticos; o efeito da reconstrução ainda não foi medido. As fórmulas históricas da Mother Track não devem ser confundidas com a implementação da versão v1.20.0 efetivamente consumida.

O valor “30” precisa sempre de unidade: percentil de que distribuição, em que versão, com que janela e normalização? Não é 30% de células acessíveis. Comparações de altura entre tracks com escalas diferentes não estabelecem diferenças biológicas.

## 5. Auditoria dos filtros e estatísticas que falta fechar

Ahmed et al. 2026 é uma revisão de critérios de safe harbor; não é um protocolo canino validado que possa ser transcrito sem decisões adicionais. [Fonte](https://pmc.ncbi.nlm.nih.gov/articles/PMC12785581/).

| Camada | Situação observada | Condição de fecho proposta |
|---|---|---|
| Universo inicial | 461 intervalos convergentes de 50–75 kb | Justificar ou ampliar de forma predefinida; documentar regiões nunca avaliadas |
| Distância a genes | Agora considera 5′ e strand; avalia a janela inteira | Fixar gene versus transcrito, extremidades, limites inclusivos e sítio local; testes sintéticos |
| Genes de risco e miRNA | Código usa raios de 300 kb | Conferir fonte, identidade/ortologia, lista completa e convenção de distância |
| RNA não codificante | Há veto de sobreposição | Conferir cobertura e versão da anotação; ausência não demonstra neutralidade |
| ATAC | p55 mais veto de pico e componente de frequência | Justificar cada papel; métricas locais por cão; medir exclusões atribuíveis |
| Repetições | Dados incompletos; corte atual >50% | Referência comparável em todos os candidatos; avaliar também inserção local |
| Conservação | Piloto e valores ausentes; corte atual phyloP >6,5 | Modelo e cobertura defensáveis, qualidade do alinhamento/liftover; limiar justificado |
| Reguladores externos | BED de três colunas sem tipos/evidência completos | Substituir/complementar por anotação rastreável; não dispensar o risco por falta de metadados |
| TAD | Contexto de um cão, célula não especificada | Justificar uso como veto; avaliar incerteza e sensibilidade; não afirmar equivalência a T cells |
| Mappability | Coluna de evidência existe | Confirmar método, comprimento de reads, assembly e qualidade local |
| RRBS | Média corrigida, cobertura reduzida | Número de CpGs, profundidade, cães elegíveis e definição da estimativa explicitados |
| Variabilidade | Ausência/baixo n pode enviesar estabilidade | Demonstrar suporte mínimo e exportar n efetivo; não chamar estabilidade a n=1 |
| Score | Ranking agregado heurístico | Riscos obrigatórios não compensáveis; separar mérito e incerteza; justificar pesos |

Os cortes numéricos acima descrevem a implementação, não recomendações já validadas para CAR-T. Falta uma matriz completa ligando passagem da literatura, interpretação, código, input e teste. O registo de exclusões deve mostrar interseções entre filtros; somar contagens individuais sobrestima o número de regiões excluídas.

A máscara soft da FASTA foi calculada nos 461, mas diverge até 16,92 pontos percentuais do RepeatMasker existente nos 43 comparáveis. Não foi usada para preencher silenciosamente os restantes. O modelo neutro phyloP de uma região de 100 kb continua piloto; um valor disponível não resolve a qualidade da estimativa.

## 6. Protocolo computacional independente proposto

Sequência de trabalho, ainda não executada integralmente:

1. Congelar objetivo CAR-T, dados incluídos, assembly e manifesto por amostra. Separar unidade biológica (cão) de biblioteca/run; documentar exclusões, cobertura e metadados disponíveis.
2. Reconstruir os inputs ATAC afetados a partir de ficheiros íntegros, verificando completude por cão e cromossoma. Reportar complexidade, duplicação, distribuição de fragmentos, TSS enrichment, FRiP e reprodutibilidade. Os indicadores seguem boas práticas de [ENCODE](https://www.encodeproject.org/atac-seq/); não importar limiares humanos sem justificar adequação à referência canina.
3. Produzir análise por cão e um agregado com definição explícita. Avaliar influência do peso QC e do gate local. Um peso que depende do próprio sinal pode alterar a população que contribui para cada bin; quantificar o efeito em vez de chamar automaticamente ao agregado uma média populacional sem viés.
4. Definir previamente a escala local de procura de inserção e a vizinhança usada para risco. Manter coordenadas BED 0-based half-open internamente e conversões explícitas para visualização.
5. Completar anotação comparável em todas as regiões que possam chegar à shortlist. Aplicar os critérios de segurança com rastreabilidade; classificar desconhecido quando a evidência obrigatória não existe.
6. Quantificar acessibilidade local, suporte entre cães, metilação observada e incerteza. Considerar a composição PBMC e covariáveis disponíveis; análises ajustadas não transformam bulk PBMC em ATAC de T cells.
7. Avaliar estabilidade por reamostragem de cães e exclusão de indivíduos influentes. Para estudar generalização, separar dados antes de escolher limiares/pesos. Comparar com controlos genómicos equivalentes em tamanho, mappability, cobertura e outras características técnicas relevantes.
8. Congelar ranking e critérios de avanço. Se nenhum candidato cumprir, aceitar resultado vazio e executar uma expansão de busca previamente definida. Não otimizar a regra para obter ANO2 ou qualquer outro nome.
9. Executar o workflow de ponta a ponta num ambiente fixado, com hashes, logs, testes de fronteira/strand/ausência e relatório automático. Uma execução repetida deve reproduzir os outputs dentro das tolerâncias documentadas.

Saída por candidato: assembly, intervalo regional, subintervalo local avaliado, genes/transcritos próximos, distâncias, flags e respetivas fontes, qualidade de cada anotação, ATAC por cão e agregado, suporte efetivo, RRBS por CpG/bin e cobertura, estabilidade estatística, limitações e decisão justificada. “Passa verificações registadas” deve continuar separado de “candidato elegível para validação CAR-T”.

## 7. Evidência CAR-T e condições de avanço experimental

Falta uma ponte específica para células T caninas: evidência no estado de interesse antes e após ativação/manipulação, quando disponível. Dados de sangue permitem priorizar; não demonstram manutenção do estado cromatínico durante fabrico, expansão ou atividade de CAR-T. Dados públicos de T cells caninas ainda não foram inventariados nesta análise; a sua existência e adequação não são presumidas.

A avaliação funcional deverá pré-especificar expressão adequada e persistente, caracterização da integração e número de cópias, alterações estruturais relevantes, efeitos locais e transcriptómicos, viabilidade, expansão e função das células T/CAR-T. Incluir controlos que permitam separar efeito do locus, da cassete e da edição. Expressão máxima não é automaticamente o objetivo: o nível deve ser compatível com o produto CAR-T pretendido. Número de amostras, duração e margens de aceitação dependem de dados piloto e desenho estatístico; não são números que se possam escolher honestamente a partir deste scoring.

Esta é uma especificação de evidência, não um SOP experimental executável nem um desenho de gRNA. Exemplos de validação funcional de safe harbors humanos ajudam a definir dimensões a medir, mas não validam coordenadas caninas por homologia: [Aznauryan et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9017210/).

## 8. Critérios de fecho, prioridades e registo de alterações

**Bloqueadores computacionais imediatos:** integridade/reconstrução ATAC; resolução regional versus local; completude e proveniência das anotações; justificação do veto ATAC; suporte RRBS/variabilidade; validação estatística de seleção; workflow integral reproduzível. Reordenar pesos antes de resolver estes pontos não fecha o protocolo.

**Fecho computacional:** regras congeladas e rastreáveis, todos os inputs obrigatórios qualificados, falhas e desconhecidos explícitos, ausência de bug conhecido com impacto não quantificado, testes pertinentes e execução completa reproduzida. Os resultados podem ser candidatos ou ausência de candidatos.

**Fecho biológico:** evidência funcional suficiente no contexto de CAR-T definido, com critérios de aceitação e limitações estabelecidos. Não é substituído por percentis, overlap entre listas ou pelo facto de um candidato ser o único sobrevivente.

Estado desta atualização: análise documental/código/manifest/resultados e literatura; não foi feito novo scoring, reconstrução ATAC, expansão dos 461 nem ensaio experimental. Preservam-se os outputs v1.21.1. A reconstrução CpG já concluída é anotação suplementar e não introduziu um novo veto nesta versão.

Ficheiros de suporte na pasta v1.21.1: `RESULTS.md`, `manifest.json`, `summary.json`, `comparison.tsv`, `candidates_scored_v1211.tsv`, `provisional_candidates.tsv`, `candidates_passing_recorded_checks.bed`, `cpg_validation.json`, `output_hashes.json`, `code/`, `test_regression.py`. Documentação histórica consultada: `04_tracks_processadas/ROS_Cfam_1.0/METHODS.md`, `05_SHIP/V3_PROGRESS_NOTES.md` e MOC Obsidian. Esses documentos históricos não são todos descrições atuais da v1.21.1.

Próxima decisão prática: concluir os bloqueadores computacionais sem depender de qualquer pessoa externa; manter o objetivo CAR-T já confirmado e distinguir claramente tarefas que exigem novos dados celulares de tarefas resolúveis com os dados locais.
