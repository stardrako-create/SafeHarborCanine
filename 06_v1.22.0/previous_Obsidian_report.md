# Safe Harbor canino para CAR-T — relatório v1.21.9

Atualizado em 14 de setembro de 2026. Esta versão acrescenta análises de conservação, identidade celular exploratória e variantes estruturais à reconstrução ATAC verificada da v1.21.8. **Há três regiões exploratórias; nenhum safe harbor canino para CAR-T está validado.**

## Conclusão para a reunião

Os problemas de cálculo identificados anteriormente foram corrigidos e a reconstrução/rescoring v1.21.8 foi verificada. Isso não resolve a questão biológica principal: acessibilidade e neutralidade funcional em células T caninas. O sangue apresenta janelas localmente acessíveis, mas o multiome tumoral disponível não fornece suporte positivo robusto nas três janelas, após uma análise de RNA mais ampla que a seleção por poucos marcadores.

Não recomendo apresentar um vencedor nem ajustar filtros para produzir ATAC mais alto. Recomendo apresentar uma shortlist auditada, distinguindo acessibilidade no sangue, evidência em T e segurança após integração. A próxima prioridade que pode alterar a decisão é evidência em T caninas relevantes, além da resolução local das anotações regulatórias.

## 1. Estado dos candidatos

O scorer v1.21.8 avaliou 461 regiões: 457 excluídas, uma insuficiente e três passes legados. A v1.21.9 não refez o scoring das 461 regiões com o novo modelo de conservação. Os scores e percentis seguintes são da v1.21.8; a nova evidência é acrescentada separadamente.

| Região | Score legado | Percentil ATAC regional | Percentil ATAC local de 1 kb | EpiC promotor/enhancer na janela de 1 kb |
|---|---:|---:|---:|---:|
| ANO2/NTF3 | 0,5073 | 72,70 | 93,98 | 800 bp |
| LOC119876429/LOC119872513 | 0,3984 | 57,43 | 96,03 | 0 bp |
| NPNT/TBCK | 0,5015 | 57,63 | 96,57 | 662 bp |

As três regiões completas têm sobreposições EpiC de promotores/enhancers: 9.616, 5.349 e 14.800 bp, respetivamente. Zero na janela LOC significa ausência nas anotações interrogadas, não ausência de função regulatória.

Os percentis usam fundos com comprimento correspondente. Uma região de dezenas de kb dilui sinais focais; o percentil de uma janela de 1 kb responde a outra pergunta. Além disso, as janelas foram selecionadas exploratoriamente por máximos anteriores. Percentil 96 não significa 96% de cromatina aberta, nem acessibilidade demonstrada em T. Alturas de tracks de experiências com normalização e processamento diferentes não são diretamente comparáveis.

Coordenadas ROS, 0-based half-open:

| Região | Janela anterior de 1 kb |
|---|---|
| ANO2/NTF3 | NC_051831.1:39728324–39729324 |
| LOC119876429/LOC119872513 | NC_051811.1:17235706–17236706 |
| NPNT/TBCK | NC_051836.1:27057411–27058411 |

## 2. Conservação: modelo publicado aplicado

Foi obtido o modelo geral de phyloP baseado em posições repetitivas ancestrais, disponibilizado no [servidor oficial UCSC/Zoonomia](https://cgl.gi.ucsc.edu/data/cactus/241-mammals-human-phylop-models/). A origem do treino está descrita no README guardado. Os 241 nomes de folhas coincidem com os do modelo piloto. Reexecutámos phyloP LRT/CONACC nas três regiões completas e projetámos os resultados para as 2.091 janelas locais através dos blocos de cadeia guardados.

| Região | Máximo da média móvel de 50 bp, região: piloto → publicado | Mesma estatística na janela de 1 kb: piloto → publicado |
|---|---:|---:|
| ANO2/NTF3 | 5,10122 → 4,37690 | 1,99976 → 1,47496 |
| LOC119876429/LOC119872513 | 4,01352 → 3,22948 | 0,97258 → 0,62494 |
| NPNT/TBCK | 4,22592 → 3,26356 | 1,56628 → 1,09474 |

Cobertura observada regional: 99,84%, 99,89% e 99,93%, respetivamente; nas três janelas anteriores, 100%. As médias móveis exigem 50 bases observadas consecutivas, sem atravessar gaps. Os resultados publicados e piloto têm os mesmos conjuntos de posições WIG.

Esta análise substitui uma calibração local frágil por uma comparação com um modelo publicado. Usa o mesmo alinhamento subjacente: não é replicação biológica independente. O limiar legado de 6,5 não foi validado para esta calibração; não foi reproduzida uma análise de FDR genómica. Valores inferiores não demonstram neutralidade funcional. A reprodução integral dos critérios do paper continua distinta desta análise de sensibilidade.

## 3. Contexto T: clustering RNA e fragmentos locais

Dados públicos [GSE244116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244116): um osteossarcoma canino primário, um dador. Após os limites QC usados anteriormente, 5.849 núcleos. A análise RNA independente usou normalização por biblioteca para 10.000, log1p, 2.000 genes variáveis por dispersão em 20 grupos de expressão média, escala com clipping ±10, PCA de 30 componentes, grafo de 20 vizinhos e Louvain em resoluções 0,5/1,0 com sementes 0/42.

Isto não reproduz o WNN multimodal dos autores nem o VST do Seurat. A anotação exploratória exigiu clusters com ≥20 núcleos, ≥30% com pelo menos um CD3D/E/G, ≥30% com CD247/LCK/CD2 e <20% com dois ou mais marcadores mieloides ou B. A acessibilidade dos candidatos não participou no clustering ou na anotação.

As quatro análises identificaram grupos enriquecidos em T de 155, 155, 157 e 157 núcleos. A interseção tem 155; a união, 157. A estabilidade desta população não implica estabilidade de todos os clusters nem identidade celular validada.

| Janela | Fragmentos com início local, interseção de 155 | União de 157 |
|---|---:|---:|
| ANO2/NTF3 | 0 | 0 |
| LOC119876429/LOC119872513 | 0 | 0 |
| NPNT/TBCK | 0 | 0 |
| Controlo CD3D | 4 | 4 |
| Controlo CD3E | 4 | 4 |
| Controlo CD247 | 7 | 7 |
| Controlo LYZ | 0 | 0 |

Os controlos são os intervalos guardados na v1.21.7; CD3D/E têm verificação adicional de TSS, CD247/LYZ são controlos centrados na feature. O quinto campo do ficheiro de fragmentos representa suporte PCR, não novas moléculas independentes. A contagem principal é por início no intervalo; não deve ser apresentada como contagem completa de todos os fragmentos sobrepostos.

Os sinais mínimos anteriores obtidos por seleção de marcadores não se mantêm nesta definição mais ampla por RNA. **Não há suporte positivo robusto para estas janelas neste grupo T-enriquecido. Zero não prova cromatina fechada:** o ensaio é esparso, a amostra é tumoral e única, e faltam deteção de doublets, QC ATAC completo e referência independente de identidade. Não são 155 dadores. Não se devem extrapolar estas contagens para T saudáveis, ativadas ou CAR-T.

## 4. NPNT: revisão da grande chamada de deleção

A chamada catalogada abrange aproximadamente 28,9 Mb no chr32 de UU_Cfam_GSD_1.0. Foram interrogados SNPs PASS em dez intervalos de 2 kb distribuídos no interior e quatro no exterior, no [VCF público Dog10K](https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/AutoAndXPAR.SNPs.vqsr99.vcf.gz). Considerámos heterozigotia com DP≥8 e GQ≥20. São 299 sítios interiores interrogados por amostra, não cobertura contínua de toda a extensão.

| Portador catalogado com FT PASS na SV | SNPs heterozigóticos de alta qualidade no interior | Intervalos interiores com esses SNPs | Mediana DP interior/exterior |
|---|---:|---:|---:|
| GBGV000001 | 10 | 5 | 19/19 |
| PBGV000010 | 22 | 7 | 19/18 |
| PBGV000012 | 26 | 4 | 19/21 |

O quarto portador, BFDB000001, tem FT FAIL3 na SV; apresentou 21 SNPs heterozigóticos de alta qualidade em cinco intervalos. Os dados de comparação por prefixo de raça também estão guardados; duas amostras têm raça desconhecida nos metadados, pelo que não são controlos de raça certificados.

A heterozigotia distribuída põe em causa uma deleção heterozigótica simples e constitutiva em toda a extensão. A profundidade nos sítios interrogados também não apresenta uma redução aproximada para metade, mas é uma medida esparsa condicionada à presença de variantes, não uma estimativa imparcial de número de cópias.

O estado correto passa a ser **chamada SV catalogada, inconsistente com perda simples de uma cópia em toda a extensão, ainda por resolver ao nível dos reads**. Não é uma deleção confirmada nem um artefacto provado. Não elimina rearranjos complexos, mosaicismo ou problemas de mapeamento. Não autoriza automaticamente NPNT.

## 5. Alternativas locais e interpretação

As 64 janelas da frente exploratória anterior foram preservadas e anotadas com o novo modelo: 22 ANO2, 12 LOC e 30 NPNT. São subintervalos frequentemente sobrepostos, não 64 loci independentes. O conjunto usa compromissos entre ATAC, cobertura RRBS, repetição e metilação observada; não foi criada uma ponderação ou um limiar de conservação para forçar um vencedor.

A tabela integrada cobre 2.091 janelas. Os zeros do multiome dizem respeito apenas às três janelas anteriores interrogadas; não se estendem às 64 alternativas. Muitas alternativas ainda não têm a revisão local de variantes, mapeamento e acessibilidade T feita para aquelas três.

LOC mantém interesse para investigação local pela ausência de sobreposição regulatória interrogada na janela anterior e menor conservação local entre as três. Isso não supera a falta de suporte T. ANO2 tem melhor percentil regional, mas a janela anterior conflita com anotação regulatória. NPNT tem o maior percentil local, mas limitações de RRBS e conflito regulatório, além da SV não resolvida. Nenhum destes argumentos permite eleger um safe harbor.

## 6. O que falta para uma decisão defensável

1. **Evidência na célula relevante:** testar/reanalisar acessibilidade em T caninas saudáveis e ativadas, com replicação entre dadores e QC declarado. O multiome atual dá contexto exploratório, não validação CAR-T. No material disponível, completar QC/identidade pode refinar a conclusão, mas não criar replicação.
2. **Auditoria local das alternativas:** antes de selecionar uma janela, aplicar-lhe as mesmas verificações de coordenadas, variantes, mapeamento, conservação e anotações regulatórias. Fixar critérios antes de comparar desempenho em dados independentes para limitar seleção por máximos.
3. **Conservação e filtros reproduzíveis:** distinguir reprodução exata do paper de alterações metodológicas justificadas; validar a decisão estatística do modelo publicado antes de o usar como veto ou refazer o scoring global. As 461 regiões não foram reavaliadas integralmente nesta versão.
4. **Sequência individual e SV:** resolver a chamada NPNT com evidência de alinhamentos/número de cópias e verificar a sequência dos animais relevantes. Catálogos populacionais e a proxy de mapeabilidade não substituem genotipagem individual ou avaliação de especificidade de edição.
5. **Validação funcional após integração:** demonstrar expressão sustentada, integridade genómica e ausência de perturbação relevante da expressão vizinha, crescimento e função T. Esses resultados experimentais faltam; nenhum ajuste computacional os substitui.

## 7. Verificação e rastreabilidade

Passaram as verificações de integração: 2.091 janelas únicas; cálculo independente por somas das médias móveis das três janelas; correspondência exata dos 5.849 barcodes com QC; recálculo das 42 contagens de fragmentos a partir dos registos em cache, sem duplicados de coordenadas/barcode; recálculo dos totais de heterozigotia e critérios DP/GQ; preservação das 64 alternativas. Estas verificações são de aritmética e associação de dados, não certificação dos experimentos de origem. Mantêm-se as condições de completude da consulta de fragmentos da v1.21.7.

Ficheiros nesta pasta: `conservation_model_comparison.tsv`, `conservation_review.json`, `conservation_execution.json`, `model_compatibility.json`, `RNA_cluster_summary.json`, `RNA_cluster_profiles.tsv`, `RNA_cluster_marker_profiles.tsv`, `RNA_cluster_assignments.tsv`, `RNA_cluster_local_fragment_counts.tsv`, `NPNT_span_genotypes.tsv`, `NPNT_span_review.json`, `local_windows_integrated_v1219.tsv`, `local_tradeoff_frontier_annotated_v1219.tsv` e `integration_validation.json`. Scripts e hashes preservados. O checkpoint v1.21.8 permanece como histórico da reconstrução, testes e rescoring.
