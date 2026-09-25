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
