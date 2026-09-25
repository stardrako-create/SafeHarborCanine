# Safe harbor canino para CAR-T — relatório de situação para reunião

**Data:** 2026-09-14 14:19 (hora local). **Estado:** resultados auditados até v1.21.7; reconstrução ATAC v1.21.8 concluída, ainda sem validação final e novo scoring. Este documento substitui o adiamento anterior do relatório, por pedido expresso do utilizador.

## Mensagem principal

O projeto produziu uma shortlist de regiões candidatas e melhorou substancialmente a consistência dos cálculos. **Ainda não identificámos um safe harbor canino validado para CAR-T.** As três regiões que passaram as verificações do scorer anterior não estão aprovadas perante toda a evidência acrescentada depois. Duas janelas locais têm anotações regulatórias independentes; a terceira continua com evidência funcional insuficiente. O conjunto público com marcadores T fornece poucas observações e não resolve esta limitação.

A formulação defensável na reunião é: **«Temos três regiões prioritárias para investigação, com limitações e riscos locais identificados. Corrigimos erros de comparação estatística e tratamento de dados em falta. Estamos a concluir a reavaliação computacional; a acessibilidade e neutralidade em células T caninas continuam por demonstrar.»**

## 1. O que foi analisado

O contexto-alvo é CAR-T canina. A evidência epigenómica principal deriva de sangue/PBMC, incluindo os 71 cães do estudo de base; a integração ATAC foi alargada a 76 amostras com cinco adicionais. Sangue misto não é equivalente a células T purificadas, ativadas ou geneticamente modificadas.

A última comparação completa anterior à reconstrução atual avaliou **461 regiões**: **457 excluídas, uma com evidência insuficiente e três que passaram as verificações então registadas**. O scoring coerente com consenso ATAC de 76 amostras reproduziu os scores e estados anteriores. Estes números referem-se ao scorer legado; não são a contagem de safe harbors validados nem o resultado de um novo scoring v1.21.8.

Depois, foram analisadas **2.091 janelas locais**, com comprimentos de 500, 1.000 e 2.000 bp, nas três regiões. Foram explorados ATAC, cobertura RRBS, metilação observada, repetições, conservação e anotações regulatórias. As três janelas de 1 kb discutidas abaixo foram escolhidas anteriormente por sinal ATAC local elevado; essa seleção exploratória favorece máximos e não constitui confirmação independente.

## 2. Resultados locais que podemos apresentar

Coordenadas na montagem **ROS_Cfam_1.0**, 0-based e fim exclusivo. Os nomes correspondem aos genes que identificam a região, não a inserções propostas dentro desses genes.

| Região | Janela local de 1 kb | Percentil ATAC local anterior | Sobreposição promotor/enhancer EpiC em 11 tecidos | SNPs + nonSNPs PASS no Dog10K |
|---|---|---:|---:|---:|
| ANO2/NTF3 | NC_051831.1:39728324–39729324 | 94,1 | 800 bp | 20 + 1 |
| LOC119876429/LOC119872513 | NC_051811.1:17235706–17236706 | 96,07 | 0 bp | 13 + 2 |
| NPNT/TBCK | NC_051836.1:27057411–27058411 | 96,63 | 662 bp | 23 + 8 |

Os percentis são dos cálculos locais anteriores, contra 3.000 janelas de fundo com comprimento correspondente, e ainda não foram recalculados com os novos tracks. Não são percentagens de cromatina aberta, nem percentis de células T, nem evidência de segurança. Não devem ser comparados diretamente com alturas de picos produzidas por outros pipelines.

**ANO2/NTF3:** a janela tem anotações de enhancer em pulmão/glândula mamária. O sinal ATAC local elevado não elimina a possibilidade de função regulatória endógena. Exige resolução do conflito antes de se defender neutralidade.

**NPNT/TBCK:** há anotação de enhancer em rim e uma deleção estrutural catalogada que atravessa a janela. A evidência RRBS local é especialmente fraca: cerca de 8,9% de cobertura observada e no máximo três cães nas posições cobertas; a média de 100% de metilação refere-se apenas à pequena fração observada. Não se pode extrapolar esse valor para toda a janela ou para CAR-T.

**LOC119876429/LOC119872513:** não sobrepõe promotor/enhancer nos 11 tecidos EpiC analisados e tem menos variantes frequentes na janela consultada. É relativamente menos problemática nestas anotações, **mas não um vencedor validado**. Os tecidos não incluem T purificadas; ausência de anotação não demonstra ausência de função. Houve ainda uma alternativa de alinhamento num dos fragmentos de 100 bp no teste local de remapeamento.

## 3. Porque apareciam percentis ATAC baixos

Há três questões distintas:

1. **A comparação estatística original tinha erros.** Algumas estatísticas calculadas em janelas eram comparadas com um fundo de bins individuais. A agregação muda a distribuição. Foram corrigidas comparações de acessibilidade, variabilidade, metilação e frequência de picos para usar estatísticas e comprimentos compatíveis.
2. **A escala espacial muda a pergunta.** Uma região de dezenas de kb pode ter média modesta e conter uma janela de 1 kb relativamente acessível. Encontrar um máximo local não garante acessibilidade ampla, robustez entre animais ou neutralidade regulatória.
3. **O pipeline tem escolhas conservadoras que também desfavorecem picos fortes.** O veto a qualquer pico consenso na região e a componente que favorece baixa frequência de picos são escolhas do projeto. Há uma tensão real entre procurar acessibilidade e penalizar picos. Não basta retirar estes critérios para obter números mais altos: é necessário distinguir acessibilidade útil de elementos regulatórios endógenos a preservar.

O artigo de Ahmed et al. de 2026 é uma **revisão sobre safe harbors humanos**, não uma validação de limiares em cães. Os valores usados no projeto devem ser apresentados como adaptações/heurísticas, quando não houver validação específica. [Revisão de Ahmed et al.](https://pubmed.ncbi.nlm.nih.gov/41511364/).

## 4. O que dizem os dados associados a células T

Foi consultado o multiome público **GSE244116**, de um osteossarcoma primário de um único Doberman. A matriz contém 5.969 núcleos; 5.849 passam os limites de contagens/features usados no script dos autores. Não reproduzimos ainda a classificação WNN dos autores. Definimos grupos exploratórios apenas com RNA, sem usar ATAC dos candidatos para escolher células: 93 núcleos com pelo menos dois genes CD3 detetados; um grupo mais restrito de 58 com marcadores adicionais; e um grupo alternativo de 133. Os grupos não são todos aninhados.

As janelas foram convertidas para canFam6 por dois percursos de chains que concordam nas 3.000 bases. Os fragmentos foram consultados regionalmente no BGZF; os intervalos foram enquadrados por registos anteriores/posteriores e a descodificação foi conferida por HTSlib. O ficheiro integral não foi descarregado. As contagens principais são registos únicos de fragmentos **com início dentro da janela**, não suporte PCR somado, nem CPM.

| Janela | Todos após QC: 5.849 | Grupo de 93 | Grupo restrito de 58 | Grupo alternativo de 133 |
|---|---:|---:|---:|---:|
| ANO2/NTF3 | 25 | 0 | 0 | 1 |
| LOC119876429/LOC119872513 | 43 | 2 | 1 | 1 |
| NPNT/TBCK | 44 | 0 | 0 | 1 |

Mesmo as janelas dos TSSs de CD3D e CD3E têm apenas dois fragmentos cada no grupo restrito. **Estes dados não sustentam acessibilidade forte em T, mas os zeros também não demonstram cromatina fechada.** Um ou dois fragmentos em LOC não estabelecem uma vantagem robusta. Faltam anotação celular completa, QC específico de ATAC e replicação biológica. Núcleos do mesmo tumor não são cães independentes. [GEO GSE244116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244116).

## 5. Variantes, conservação e outras limitações

**Variantes pequenas:** foram recuperados 78 registos nas três janelas. Todos os alelos REF coincidem com o FASTA UU usado na consulta, e duas implementações independentes devolveram as mesmas contagens. Com AF alternativa ≥1%, há 9 registos PASS em ANO2, 3 em LOC e 12 em NPNT. Estes números combinam SNPs/nonSNPs e não são um veto automático. O release tem 1.987 amostras; não corresponde aos 76 cães do projeto nem ao subconjunto de 1.929 da publicação. As outras 2.088 janelas não foram avaliadas para variantes. [Dog10K SNP/nonSNP](https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/).

**SV em NPNT:** uma deleção agregada de aproximadamente 28,9 Mb, AF de cerca de 0,1064%, sobrepõe a janela. Há quatro genótipos heterozigóticos catalogados, três FT PASS, em raças basset. Nos três PASS não há heterozigotia nos 23 SNPs locais analisados, mas homozigotia pode ter outras explicações. Não confirmámos a deleção nem a sua presença nos animais do projeto. Mantém-se como alerta que requer evidência de leituras ou validação ortogonal. [Catálogo SV e critérios](https://kiddlabshare.med.umich.edu/dog10K/Manta-SV_2022-03-28/svcallset.v2.README.txt).

**Conservação:** o modelo permanece piloto. Reajustar apenas a escala da árvore com uma segunda região aleatória mudou-a para cerca de 85% da original e alterou os valores locais de phyloP. Essa região não foi filtrada como conjunto neutro de referência, e o modelo herda parâmetros do piloto. Não podemos declarar a conservação calibrada independentemente. A reprodução local diferiu até 0,034/0,033 por base em LOC/NPNT; causa completa por resolver.

**QC e estabilidade:** a análise por cão foi feita nas três janelas e regiões, sem erros nessas leituras locais. O bootstrap emparelhado dos 71 cães não separa claramente os candidatos: os intervalos das diferenças incluem zero. O ranking varia com ponderação/QC. Só 2/76 atingem FRiP ≥0,2 no cálculo disponível, mas não foi estabelecido um critério universal de aprovação. O TSS enrichment do pipeline usa uma definição que não equivale diretamente à métrica de pico ENCODE.

**Hi-C e contexto celular:** o dado disponível de um cão não resolve a arquitetura cromatínica em T/CAR-T. Distâncias e fronteiras TAD são evidência contextual, não prova de neutralidade. A anotação EpiC de 11 tecidos acrescenta conflitos regulatórios que o score anterior não refletia integralmente. [EpiC Dog](https://doi.org/10.1126/sciadv.ade3399).

## 6. Atualização de hoje: v1.21.8

A reconstrução global ATAC **terminou a execução**: 76 amostras, 376 sequências, 2.396.858.295 bp e bins de 25 bp; duração registada de aproximadamente 23,3 minutos. Quatro intervalos de 2 Mb foram recuperados a partir dos BigWigs raw convertidos pela escala CPM, em três amostras. Foi produzida uma track adicional com número de cães que passam o gate por bin. Bins sem cães com evidência ficam sem média observada; variabilidade exige pelo menos dois cães.

Foram também alterados dois comportamentos no código principal: falhas de leitura já não podem ser convertidas silenciosamente em zero, e componentes de score requeridas mas ausentes passam a constar como evidência em falta, mesmo que haja score parcial calculável.

**O que esta atualização ainda não significa:** não há novo ranking aprovado. Falta validar a recuperação raw→CPM, comparar os novos tracks com medições independentes, testar especificamente as alterações novas, recalcular scores e rever diferenças. Os 11 testes de regressão existentes passaram hoje; esses testes não cobrem integralmente o novo construtor, as recuperações ou a nova regra de componentes ausentes. Portanto não declarar «pipeline totalmente fechado».

## 7. O que falta, por ordem de decisão

| Prioridade | Trabalho em falta | Resultado necessário |
|---|---|---|
| Imediata, computacional | Validar v1.21.8, recuperações e alterações; recalcular scoring e rever diferenças | Uma versão de resultados consistente, rastreável e revista |
| Imediata, seleção | Resolver conflitos regulatórios e avaliar janelas alternativas sem otimizar apenas o percentil ATAC | Shortlist justificada para investigação, com riscos e missing explícitos |
| Alta, contexto celular | Anotação/QC completos do multiome e evidência em T caninas saudáveis/ativadas de vários animais | Suporte reprodutível de acessibilidade no tipo celular relevante |
| Alta, genómica | Conservação calibrada, SV NPNT e variantes da sequência individual | Reduzir incerteza estrutural e de sequência |
| Indispensável, experimental | Demonstrar expressão estável e ausência de perturbação relevante após integração, no contexto CAR-T | Validação funcional de safe harbor; não é substituível por score |

Ainda não foi escolhido nem validado um ponto de inserção individual, e a análise local de remapeamento não equivale a análise completa de off-target. A literatura funcional em T humanas ilustra a necessidade de testar estabilidade e função após integração; os resultados não se transferem automaticamente para cães. [Odak et al., Blood](https://pmc.ncbi.nlm.nih.gov/articles/PMC10273162/).

## 8. Afirmações a evitar na reunião

- «Temos três safe harbors caninos» — temos três regiões da shortlist legada, com conflitos e lacunas posteriores.
- «Percentil 96 significa altamente acessível em CAR-T» — refere-se a uma estatística local de outra população de dados, escolhida exploratoriamente.
- «Ausência de pico, de anotação ou de cobertura prova segurança» — são situações distintas; missing não é zero biológico.
- «Um overlap é um sinal real» — para 26 e 19 elementos num universo de 461, a expectativa simples por acaso é 26×19/461 ≈1,07. Um overlap isolado não estabelece enriquecimento; o nulo depende também do desenho da seleção.
- «A reconstrução terminada valida os resultados finais» — a execução terminou, mas validação e novo scoring estão pendentes.

## 9. Rastreabilidade

Dados e verificações estão nas pastas 06_v1.21.4 a 06_v1.21.8. Destacam-se: v1.21.5/cohort_and_support.tsv e local_windows_regulatory_rechecked.tsv; v1.21.6/variant_validation.json, NPNT_SV_carriers_review.json e conservation_sensitivity_summary.json; v1.21.7/local_T_marker_fragment_counts.tsv, mapping_provenance.json e fragment_quantification_validation.json; v1.21.8/ATAC/build_manifest.json e recovery_events.json. Este relatório é um retrato do estado da reunião, não uma certificação do pipeline nem validação biológica.
