# Safe Harbor CAR-T — v1.22.4: contexto de expressão em CAR-T caninas

16 de setembro de 2026. **A matriz de RNA-seq CAR-T canina GSE247355 foi analisada descritivamente e associada ao contexto génico dos candidatos. Há genes vizinhos expressos nas células T/CAR-T; a ausência de alguns símbolos na tabela processada não demonstra silêncio. Não foram alterados os filtros, scores ou nomeações.**

## Desenho e auditoria de amostras

O [GEO GSE247355](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247355), associado a [Cao et al., 2024](https://pubmed.ncbi.nlm.nih.gov/38554158/), contém nove colunas: células T não transduzidas, CAR B7-H3 e CAR B7-H3/CXCR2, cada condição com os rótulos de réplica biológica B, E e M. As nove colunas foram reconciliadas individualmente com GSM7887650–GSM7887658 através do ficheiro SOFT oficial. Títulos, tratamento e rótulo de réplica concordam com a matriz. Não se contou isto como nove dadores independentes, nem se confirmou identidade individual para além dos rótulos depositados.

Há quatro conflitos explícitos nos campos source_name/tissue: GSM7887652, GSM7887653, GSM7887654 e GSM7887655 referem outra condição ou réplica nesses campos. Duas amostras têm esses campos específicos concordantes; três têm origem genérica T cell e tecido não preenchido. A tabela preserva todos os valores e marca cada situação. Por exemplo, GSM7887652 tem título B_T_cell_1, réplica B e tratamento Untransduced, mas source_name/tissue E_B7H3_CAR_T. Usámos os títulos/characteristics concordantes para descrever a matriz; a identidade experimental subjacente não foi revalidada por genótipos.

O protocolo depositado indica análise de culturas nos dias 10–14 e processamento contra CanFam3.1. O controlo é T não transduzida no contexto deste estudo, não um substituto automático de células T quiescentes recém-isoladas. A biblioteca é RNA; não fornece ATAC nas janelas candidatas.

## Resultado de expressão

Valores abaixo são medianas CPM das três réplicas rotuladas por condição, calculadas da matriz depositada; não são testes de expressão diferencial.

| Gene | T não transduzida | CAR B7-H3 | CAR B7-H3/CXCR2 |
|---|---:|---:|---:|
| NPNT | 0,0731 | 0,0357 | 0,1499 |
| TBCK | 17,0276 | 15,2804 | 17,3845 |
| TSEN15 | 8,0881 | 8,9919 | 8,0444 |
| C7H1orf21 | 6,6905 | 6,8129 | 7,9989 |
| ARPC5 | 148,946 | 134,206 | 123,122 |
| TET2 | 93,2744 | 92,5399 | 87,3920 |

TBCK, TSEN15 e C7H1orf21 estão detetados nas nove colunas. Isto impede descrever o contexto dos candidatos como uma vizinhança globalmente inativa em CAR-T. Não significa que a inserção proposta perturbe esses genes nem estabelece contacto regulatório à distância.

NPNT tem valores muito baixos nesta matriz, de zero a 0,232385 CPM nas nove amostras. Não foram convertidos em prova de ausência de função, permissividade à integração ou abertura do intervalo intergénico. ARPC5 e TET2 são genes de contexto/risco mais distantes; a expressão apresentada não demonstra que sejam regulados pelas janelas.

ANO2, NTF3, LOC119876429 e LOC119872513 não aparecem sob esses símbolos na tabela processada. Também faltam símbolos de marcadores como CD247 e CD8B. Isso reforça que a tabela de 14.385 genes não permite transformar ausência de linha em ausência de expressão: podem intervir anotação, símbolos, filtros ou processamento. Resolver cada ausência requer rastrear a anotação utilizada e, se necessário, reanalisar contagens/reads.

## Âmbito e validação

Selecionámos dez símbolos génicos distintos mais próximos de cada uma de três âncoras declaradas, além dos genes candidatos e marcadores já definidos, sem escolher genes pelo valor de expressão. As âncoras são ANO2 39728324–39729324, LOC 17255706–17256706 e NPNT 27041161–27042161 nos respetivos cromossomas ROS. São referências espaciais de contexto, não novas nomeações. A lista final tem 40 genes, 25 presentes e 15 ausentes por correspondência exata de símbolo. Coordenadas, cadeia e distâncias constam da tabela.

Foram produzidos intervalos mínimo–máximo e medianas por condição, e razões log2 dentro de cada rótulo B/E/M quando os dois valores são positivos. Não foi adicionado pseudocount; zeros e genes ausentes têm estados distintos. As razões são apenas descritivas. Não se fizeram testes DESeq2 com CPM, nem p-values, inferência de equivalência ou classificação de genes como estáveis. A inexistência de grandes diferenças visuais não prova estabilidade funcional.

Uma segunda leitura da matriz confirmou todas as medianas e extremos apresentados. Foram verificados nove títulos únicos associados às nove colunas e a concordância dos tratamentos/réplicas; quatro conflitos source/tissue foram preservados em vez de ignorados. A verificação é de associação e aritmética, não de identidade biológica independente.

## Consequências para o projeto

1. O dataset oferece contexto diretamente em CAR-T caninas e deve acompanhar a shortlist. O multiome tumoral deixa de ser a única referência celular discutida, mas continua a ser a fonte usada aqui para acessibilidade local exploratória.
2. Não há suporte para chamar estas regiões biologicamente neutras pela baixa expressão de um dos genes que lhes dá nome. Há outros genes próximos expressos, e os nomes históricos das regiões não substituem a anotação atual.
3. Permanecem abertos o ATAC em T/CAR-T relevantes, identidade/QC completo do multiome e neutralidade após integração. Este RNA-seq não mediu a consequência de inserir o CAR nos nossos loci: não é um ensaio de segurança desses locais.
4. A próxima lacuna documental é resolver símbolos/genes ausentes e procurar contagens não normalizadas/anotação exata; a próxima lacuna de identidade celular pode beneficiar dos atlas caninos saudáveis já identificados. A discordância dos metadados deve ser resolvida antes de inferências entre condições que dependam da identidade amostral.

Ficheiros: `sample_metadata_audit.tsv`, `GEO_sample_metadata.json`, `GSE247355_family.soft.gz`, `local_gene_context.tsv`, `CAR_T_gene_expression_context.tsv`, `CAR_T_paired_descriptive_ratios.tsv`, `RNA_context_validation.json` e scripts. O scoring global continua v1.21.8. **Nenhum safe harbor CAR-T validado.**


---

# Histórico anterior até v1.22.3

A análise descritiva do RNA CAR-T, anteriormente pendente, está apresentada acima; análise inferencial e identidade individual não foram concluídas.

# Safe Harbor CAR-T — v1.22.3: profundidade do multiome e dados públicos CAR-T caninos

15 de setembro de 2026. **Há dados públicos diretamente de CAR-T caninas. A referência RNA GSE247355 foi agora identificada, descarregada e verificada. Além disso, a análise da profundidade mostra que os zeros do multiome tumoral têm pouca força para argumentar contra acessibilidade.** Nenhum destes resultados demonstra um safe harbor.

## Dados CAR-T caninos encontrados

[GSE247355](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247355), associado ao [estudo de Cao e colaboradores, 2024](https://pubmed.ncbi.nlm.nih.gov/38554158/), deposita RNA-seq de CAR-T B7-H3, CAR-T B7-H3/CXCR2 e células T controlo. Há nove amostras/colunas, organizadas em prefixos B, E e M. A identidade e condições dos dadores ainda não foram reconciliadas individualmente nesta análise; não se deve inferir automaticamente o desenho apenas dos nomes.

A matriz CPM pública foi descarregada: 14.385 linhas génicas e nove colunas numéricas, todas finitas e não negativas. Os reads brutos estão associados ao SRA. Foi criada uma tabela de presença e valores para genes próximos dos candidatos e marcadores T, como inventário inicial. Ausência de um gene nesta tabela processada não foi convertida em zero. Ainda não foi feita análise diferencial; CPM não deve ser entregue como contagens brutas a um modelo que exija counts.

Este dataset é diretamente útil para contexto de expressão em células CAR-T caninas. Não é ATAC-seq e não resolve diretamente se a janela de inserção está aberta. Também não demonstra neutralidade de integração nos nossos loci. A sua ausência nas iterações anteriores era uma lacuna da pesquisa, não ausência geral de dados CAR-T na literatura.

O estudo [Unedited allogeneic iNKT cells show extended persistence in MHC-mismatched canine recipients](https://pmc.ncbi.nlm.nih.gov/articles/PMC10591065/) disponibiliza GSE229457 (scRNA-seq) e GSE229458 (Nanostring). Inclui investigação de iNKT e engenharia CAR; são recursos complementares de outra população celular. Não se deve classificar todo o material como CAR-T convencional, nem chamar RNA-seq ao painel Nanostring. As condições de cada amostra precisam de verificação antes de integração.

Para identidade celular canina independente, também foram identificados [GSE225599](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225599), com leucócitos de sete cães saudáveis e dez com osteossarcoma, e [GSE301630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301630), com seis cães saudáveis. São referências transcriptómicas, não ATAC CAR-T. Nesta pesquisa não foi confirmado um dataset público de ATAC-seq especificamente de CAR-T caninas. Isto descreve o resultado da pesquisa, não prova inexistência.

## Profundidade do grupo T no multiome usado

Recalculámos diretamente do H5 as métricas RNA/ATAC da matriz e confirmámos a concordância exata com o cache de todos os 5.969 barcodes. No conjunto QC de 5.849 núcleos, os 155 núcleos T consensuais representam 2,65% dos núcleos, mas apenas **1,07% das contagens da matriz de picos**. A mediana é 2.787 contagens por núcleo T, contra 6.064 nos restantes núcleos QC excluindo os 157 da união T.

Isto confirma uma diferença de profundidade da matriz. Não é total de fragmentos genómicos, FRiP, TSS enrichment ou um certificado de qualidade. Continuam por avaliar adequadamente nucleossomas, doublets e identidade independente.

## Os zeros são informativos?

Fizemos 2.000 amostragens descritivas, com semente 42, de grupos de 155 núcleos dos restantes QC, preservando a distribuição dos T em dez estratos de profundidade ATAC da matriz. Foram excluídos os 157 da união T do conjunto de comparação. Cada amostragem escolhe sem reposição dentro dos estratos; os grupos entre amostragens podem reutilizar núcleos. A profundidade agregada mediana destes grupos foi 559.883,5 contagens, próxima das 550.644 dos T.

Nas três janelas originais, os grupos de comparação com profundidade semelhante tiveram zero inícios de fragmentos em **63,45% (ANO2), 44,40% (LOC) e 41,45% (NPNT)** das amostragens. Portanto, zero nos T não é, por si só, um resultado surpreendente nesta escala de deteção. Isto reforça a conclusão já cautelosa: falta suporte positivo robusto, mas não há demonstração de cromatina fechada.

Nos controlos CD3D/CD3E/CD247, observaram-se 4/4/7 inícios nos T, enquanto as medianas dos grupos de comparação foram zero. É coerente com a anotação exploratória T, mas não valida pureza, identidade individual ou acessibilidade de todos os candidatos.

As distribuições de reamostragem são comparações condicionais dentro de um tumor, não intervalos de confiança entre dadores ou testes de significância biológica. O grupo de comparação inclui populações heterogéneas. A profundidade por picos é uma proxy e não corrige todos os vieses. A tabela cobre as 37 alternativas com mapeamento integral/contínuo e os sete intervalos originais/controlos; não estende a comparação integral às 27 alternativas parcialmente mapeadas.

## Prioridade atualizada

A próxima análise deve incorporar a expressão em CAR-T caninas de GSE247355 e usar referências caninas saudáveis para melhorar a identidade celular do multiome. O multiome tumoral permanece evidência auxiliar de acessibilidade esparsa. Nem os seus zeros nem a expressão de genes vizinhos devem virar um veto ou um passe automático.

Ficheiros: `depth_group_summary.tsv`, `depth_matched_local_counts.tsv`, `depth_review_summary.json`, matriz CPM GSE247355, `CAR_T_source_inventory.json`, `CAR_T_context_gene_inventory.tsv`. A pesquisa de fontes e o inventário não equivalem a uma análise completa dos novos estudos. A shortlist e o scorer global não foram alterados.


---

# Histórico anterior até v1.22.2

As fontes CAR-T e a avaliação da força dos zeros foram atualizadas acima.

# Safe Harbor CAR-T — v1.22.2: correspondência LOC e detalhe de mapeabilidade NPNT

15 de setembro de 2026. **A correspondência central da janela LOC pendente ficou sustentada por alinhamento recíproco, permitindo consulta suplementar de variantes. Os flancos apresentam complexidade adicional. A investigação dos alinhamentos NPNT esclareceu que são alternativas imperfeitas, não cópias exatas equivalentes. Nenhuma janela foi promovida a safe harbor.**

## LOC: o que foi resolvido

Janela ROS NC_051811.1:17234956–17235956, 0-based half-open, correspondente ao intervalo UU NC_049228.1:17382280–17383266 (986 bp). O alinhamento original tem CIGAR `452M16I183M2D349M`, MAPQ60 e NM23. A decomposição conferida contra os FASTAs é de cinco substituições, 16 bases ROS sem par UU e duas bases UU sem par ROS.

O alinhamento da sequência UU de volta a ROS regressa exatamente às fronteiras da janela de 1 kb, MAPQ60, sem alternativa reportada, CIGAR `452M16D183M2I349M`. Os **984 pares de coordenadas são idênticos nos dois sentidos**, sem pares discordantes. Isto sustenta a correspondência desta janela central. Não demonstra igualdade entre as duas sequências nem genótipo de um animal experimental.

A regra operacional anterior NM≤10 impediu a consulta automática, apesar desta correspondência central. Foi mantida no resultado histórico, mas deixou de ser tratada como impedimento absoluto à obtenção de informação: fizemos uma consulta suplementar identificada como tal. Não foi baixado um filtro biológico nem alterado o scorer global.

## LOC: resultado suplementar de variantes

No intervalo UU correspondente, os [callsets públicos Dog10K](https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/) de 1.987 amostras devolveram:

| Callset | Registos | PASS | PASS com AF alternativa ≥1% |
|---|---:|---:|---:|
| SNP | 22 | 19 | 7 |
| nonSNP | 1 | 1 | 1 |

Todos os REF conferem com o FASTA UU. A leitura independente do texto VCF reproduz os totais, filtros e contagens de AF. Nenhuma SV sobreposta foi encontrada no catálogo SV consultado para esse intervalo; isso não prova ausência de SV.

Dos 23 registos, um SNP não tem projeção do intervalo REF contínua/disponível em ROS. Para os restantes, a tabela fornece apenas a projeção das coordenadas do REF. **Os alelos ALT não foram convertidos nem normalizados para ROS.** Mesmo um intervalo REF projetável não autoriza copiar automaticamente o alelo entre referências que diferem. As 16 bases ROS sem par UU continuam sem avaliação de SNPs por este VCF de referência UU; não são bases comprovadamente invariantes.

## LOC: os flancos não ficam validados

Também interrogámos a janela com 1 kb adicional de cada lado. Tanto no sentido ROS→UU como no inverso há alinhamentos divididos/suplementares. Os alinhamentos principais deixam 717 bases no início sem alinhar; no sentido inverso há ainda alternativas noutros contigs. Os SAMs completos foram guardados.

Assim, a reciprocidade exata da janela central **não se estende a toda a vizinhança de aproximadamente 3 kb**. A causa dos segmentos divididos pode envolver diferenças entre montagens e sequência repetitiva; não foi determinada nesta etapa. Este é um motivo para manter revisão de sequência local antes de qualquer desenho experimental. Não se inferiu uma alteração estrutural no dador.

## NPNT: natureza dos alinhamentos alternativos

Na janela NC_051836.1:27041661–27042661, os oito segmentos que falharam a proxy anterior abrangem, no conjunto, **NC_051836.1:27042361–27042586**. São quatro segmentos de 100 bp e quatro de 150 bp.

Foram registados oito alinhamentos primários exatos à origem e **31 alinhamentos alternativos reportados**. Estes últimos são parciais/imperfeitos, têm MAPQ máximo 0 e NM mínimo 0. Logo, o resultado anterior não demonstra cópias exatas igualmente boas. Demonstra que alguns segmentos têm semelhança com outros locais e falham a proxy conservadora de ausência de qualquer alternativa reportada.

A extensão indicada é o intervalo abrangido pelos segmentos testados, não a fronteira mínima de uma repetição. Retirar essas bases não foi demonstrado preservar o sinal T ou tornar a janela adequada. Não foi criada uma janela nova para evitar o aviso, nem realizada análise de especificidade de edição.

## Interpretação integrada

As 16 alternativas têm agora revisão local, incluindo consulta suplementar da única janela anteriormente não elegível, com limites de projeção explícitos. A janela LOC central deixou de ter correspondência por esclarecer, mas os seus flancos e as bases sem par permanecem pendências distintas. NPNT conserva o aviso de mapeabilidade, agora descrito sem confundir alternativas imperfeitas com duplicações exatas.

O resultado continua exploratório: o sinal T é escasso, vem de um único tumor e depende da definição celular. Ainda faltam QC/identidade completos do multiome, evidência replicada em T caninas relevantes, sequência individual e validação funcional após integração. O scoring das 461 regiões permanece v1.21.8; não foi reexecutado nesta etapa.

Ficheiros: `LOC_mapping_resolution.json`, `LOC_supplementary_lookup.json`, `LOC_supplementary_variants.tsv`, VCFs suplementares, `NPNT_failed_tile_alignments.tsv`, `NPNT_mapping_review.json`, `NPNT_alternative_summary.json`, `alternatives_integrated_v1222.tsv`, SAMs/logs de reciprocidade e flancos, `validation_summary.json`. As colunas históricas foram conservadas e a evidência nova acrescentada separadamente.


---

# Contexto histórico até v1.22.1

A atualização acima resolve a correspondência central LOC e acrescenta a consulta suplementar; pendências históricas devem ser lidas à luz dessa atualização.

# Safe Harbor CAR-T — v1.22.1: auditoria local das alternativas com sinal T

14 de setembro de 2026. **Revisão concluída de mapeabilidade e contexto anotado das 16 alternativas; variantes consultadas nas 15 que satisfazem o critério de remapeamento. Nenhuma é um safe harbor validado.** Foram incluídas as 13 janelas com sinal no grupo consensual de 155 núcleos e as três janelas ANO2 com sinal apenas na união de 157. Continuam a ser sinais de um único tumor, em janelas parcialmente sobrepostas.

## Resultado que altera a interpretação

As alternativas não são sequências invariantes nem passam uniformemente os testes de mapeabilidade. Uma janela LOC necessita de resolução adicional do remapeamento entre referências; uma janela NPNT tem alinhamentos alternativos em segmentos curtos. Nas janelas consultáveis existem variantes populacionais, incluindo variantes com frequência alélica alternativa de pelo menos 1%. Não foi criado um limiar de carga de variantes para eliminar candidatos ou promover outros.

Os três exemplos discutidos anteriormente ficam assim:

| Janela ROS, 0-based half-open | Inícios T consensuais | SNPs PASS | nonSNPs PASS | Registos PASS com pelo menos uma AF alternativa ≥1% | Proxy de mapeabilidade |
|---|---:|---:|---:|---:|---|
| LOC NC_051811.1:17255706–17256706 | 2 | 24 | 0 | 10 | 103/103 segmentos passam |
| NPNT NC_051836.1:27041161–27042161 | 2 | 16 | 2 | 9 | 103/103 passam |
| NPNT NC_051836.1:27050411–27051411 | 4 | 23 | 3 | 13 | 103/103 passam |

nonSNP é o nome da classe do callset, não uma garantia de que cada registo é uma simples inserção/deleção. As AF descrevem alelos alternativos na população do catálogo; não são probabilidades de falha nem genótipos dos animais experimentais. As contagens em janelas sobrepostas não são independentes.

LOC mantém interesse para aprofundar o contexto local, mas duas observações T e dez registos populacionais frequentes não justificam uma escolha final. NPNT não ganha prioridade só por ter quatro observações T: essa janela tinha 490 inícios no conjunto QC e apenas 3,9% de RRBS observado. Não foi demonstrado enriquecimento específico em T.

## Mapeamento entre referências e variantes

As 16 sequências ROS de 1 kb foram alinhadas à referência UU_Cfam_GSD_1.0 com BWA-MEM, incluindo saída de alinhamentos alternativos. Conservou-se a regra anterior para consulta: alinhamento integral, MAPQ≥30, ausência de alternativas reportadas, NM≤10 e alias de contig único por comprimento exato da referência. Quinze passaram. Esta é uma regra operacional de consulta, não um filtro biológico validado de safe harbor.

A janela LOC NC_051811.1:17234956–17235956 tem MAPQ60 e alinhamento integral, mas NM23, CIGAR `452M16I183M2D349M`. A diferença inclui indels entre referências; não deve ser confundida com 23 SNPs populacionais ou prova de alinhamento errado. Como excede a regra mantida, os campos de variantes ficam ausentes, não zero. Para a recuperar, é necessário justificar a correspondência ao nível dos blocos/sequências e rever a regra de consulta de forma explícita, sem a promover automaticamente.

Nas quinze janelas elegíveis consultámos os callsets públicos [SNP e nonSNP Dog10K](https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/), cada um com 1.987 amostras. Foram recuperadas 404 associações janela–registo, correspondentes a **241 registos distintos** entre callsets e coordenadas. O REF de todas as associações foi confirmado contra o FASTA UU. Os VCFs regionais preservam cabeçalhos, genótipos e FILTER originais.

Uma leitura independente do texto VCF reproduziu as contagens por sobreposição do intervalo REF, PASS e AF≥1%, sem usar os objetos pysam usados na extração. Isto valida as contas e a associação às janelas, mas não é uma segunda prova independente da completude do índice remoto.

No catálogo SV consultado, as nove alternativas NPNT elegíveis sobrepõem a mesma grande chamada de deleção anteriormente revista; não são nove SVs independentes. Nas seis janelas elegíveis ANO2/LOC não foram encontrados registos SV sobrepostos nesse catálogo. A ausência de registos não demonstra ausência de variantes. Mantém-se a interpretação NPNT da v1.21.9: chamada questionada por heterozigotia distribuída, não confirmada como perda simples de uma cópia e ainda por resolver com reads.

## Mapeabilidade local

Foram alinhados contra ROS segmentos de 100, 150 e 250 bp, com passo de 25 bp, extraídos das 16 janelas: **1.648 testes de segmentos**, incluindo segmentos repetidos entre janelas sobrepostas. O critério exige alinhamento exato à origem, MAPQ≥30 e ausência de alinhamento alternativo reportado. Passaram 1.640.

As oito falhas pertencem à janela NPNT NC_051836.1:27041661–27042661: quatro segmentos de 100 bp e quatro de 150 bp têm alternativas reportadas. Os MAPQ mínimos desses grupos são 43 e 60. Logo, MAPQ elevado isoladamente não resolve a questão. A janela deve manter um aviso de mapeabilidade local; não foi declarado todo o intervalo inutilizável.

Estes testes são uma proxy dependente do alinhador, não uma prova exaustiva de unicidade e não uma avaliação de off-target de edição. As restantes 15 janelas não apresentaram falhas sob esta proxy.

## Contexto regulatório e genes próximos

Recalculámos a sobreposição usando os blocos EpiC ROS guardados e os BEDs do projeto. Todas as 16 janelas têm zero sobreposição de promotor/enhancer EpiC, do BED regulatório externo e dos BEDs de lncRNA/smallRNA e miRNA interrogados. Isto era parcialmente esperado, porque a frente de alternativas foi pré-selecionada pela ausência de sobreposição regulatória. Não é evidência independente de neutralidade. As anotações EpiC abrangem onze tecidos, sem T purificadas.

As distâncias foram recalculadas à fronteira do corpo génico e à extremidade 5′ orientada pela cadeia, mantendo os nomes históricos das regiões separados dos genes mais próximos na anotação atual. Por exemplo, na janela LOC 17255706–17256706, o gene mais próximo nesta anotação é C7H1orf21, com intervalo de 110.509 bp até à fronteira/5′; o nome histórico LOC119876429/LOC119872513 não deve ser usado como se identificasse necessariamente os genes mais próximos nesta tabela.

As três janelas ANO2 positivas apenas na união estão a 8.427–8.927 bp da fronteira de NTF3, embora a extremidade 5′ anotada esteja a cerca de 80 kb. Distância ao corpo génico e distância ao promotor são critérios distintos; a tabela apresenta ambos. Usámos uma extremidade 5′ por registo génico, não uma análise exaustiva de todos os TSS alternativos. O catálogo de genes de risco é a lista existente do projeto, não uma avaliação exaustiva de oncogenicidade.

## Decisão e trabalho que permanece

Não há suporte para promover uma destas janelas a safe harbor. A auditoria reduziu a incerteza técnica e identificou limitações específicas, sem alterar o scorer global da v1.21.8 ou os seus passes legados.

Antes de uma escolha final faltam: resolver as correspondências de sequência problemáticas; completar identidade/QC ATAC do multiome; confirmar acessibilidade em T caninas relevantes com replicação; verificar sequência/genótipo dos animais a usar; e demonstrar estabilidade de expressão e neutralidade funcional após integração. A revisão local agora feita permite orientar essas etapas; não substitui os dados experimentais.

Ficheiros auditáveis: `alternatives_integrated_v1221.tsv`, `review_windows.tsv`, `dog10k_mapping.tsv`, `small_variant_records.tsv`, `small_variant_summary.tsv`, VCFs regionais, `SV_records.tsv`, `tile_mappability_proxy.tsv`, `tile_mappability_summary.tsv`, `regulatory_context.tsv`, `epic_state_hits.tsv`, `mapping_execution.json`, `variant_query_provenance.json`, `regulatory_provenance.json` e `integration_validation.json`. Scripts e hashes preservados.


---

# Contexto anterior: v1.22.0 e v1.21.9

A revisão local descrita acima atualiza a pendência de variantes/mapeabilidade das alternativas deste histórico.

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


---

# Contexto anterior preservado — relatório v1.21.9

A análise acima atualiza a acessibilidade das alternativas, que neste checkpoint anterior ainda estava pendente.

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
