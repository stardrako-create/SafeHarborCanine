# Auditoria adicional de fecho — em curso, 19 setembro 2026

**O painel anterior é provisório, não encerra o trabalho solicitado.** O utilizador pediu continuação e updates de cinco minutos. Este checkpoint documenta trabalho novo, não apenas reclassificação dos resultados anteriores.

## Variantes resolvidas entre referências

Confirmada igualdade integral das duas sequências de 1 kb entre ROS e UU, incluindo reverse complement para w11. Convertidos 46 registos do catálogo, incluindo FILTER não-PASS preservados; normalização com bcftools1.24 produziu48 alelos. Todos os REF conferem com ROS e o haplótipo alternativo antes/depois da normalização é idêntico. A normalização recuperou bases minúsculas do FASTA softmasked; a validação é insensível a maiúsculas, mantendo a sequência. A versão inicial do validador parou nessa diferença de caixa e foi corrigida antes de aceitar os resultados.

Ficheiros: variant_projection_provenance.tsv, projected_unnormalized.vcf, projected_normalized_ROS.vcf, normalized_variant_alleles_ROS.tsv, variant_masks_PASS_ROS.tsv. AF guardada na tabela é o máximo por registo fonte, **não a AF específica de cada alelo desdobrado**. Não copiar esse máximo para uma conclusão de frequência alélica individual. As máscaras de REF não representam toda a área de efeito de um indel. Não foram inferidos genótipos individuais.

## Especificidade de sequência em todo o assembly

Percorridos 2,396,858,295 bp em 376 contigs. Todos os segmentos consecutivos de20,25,50,100 bp nas duas janelas foram procurados exatamente nas duas orientações. Teste do algoritmo com sobreposições e fronteiras de blocos concordou com enumeração direta num exemplo controlado.

| Janela | Tamanho | Segmentos únicos | Segmentos com outras ocorrências exatas | Máximo de loci |
|---|---:|---:|---:|---:|
| w01 | 20 | 962/981 | 19 | 3 |
| w01 | 25 | 976/976 | 0 | 1 |
| w01 | 50 | 951/951 | 0 | 1 |
| w01 | 100 | 901/901 | 0 | 1 |
| w11 | 20 | 956/981 | 25 | 38 |
| w11 | 25 | 970/976 | 6 | 12 |
| w11 | 50 | 951/951 | 0 | 1 |
| w11 | 100 | 901/901 | 0 | 1 |


Isto fecha uma lacuna real da proxy antiga de100–250 bp: w01 tem19 segmentos de20 bp não únicos; w11 tem25 de20 bp e6 de25 bp. Não significa que todos sejam guias viáveis, que todos tenham PAM ou que sejam off-targets reais. Ainda falta análise de semelhanças com mismatches/bulges no contexto da nuclease/guia. Nenhuma sequência de guia foi escolhida. Contigs alternativos podem aumentar o número de ocorrências; todos os contigs do FASTA estão incluídos e os primeiros10 locais de cada padrão não único foram guardados.

## Regulação: corrigir o âmbito de “gene mais próximo”

O ficheiro histórico canine_all_genes_stranded.bed é de genes codificantes; os lncRNAs estavam noutro BED. As distâncias anteriores devem ser designadas **distâncias a genes codificantes naquele BED**, não distância a qualquer gene. GFF completo e RNA filho de gene foram agora usados para distinguir corpos e5′ por biotipo.

- w01: corpo génico mais próximo considerando todos os biotipos = LOC119872513 (17,805 bp); 5′ de transcrito codificante mais próximo = C7H1orf21 (110,509 bp).
- w11: corpo génico mais próximo considerando todos os biotipos = NPNT (27,250 bp); 5′ de transcrito codificante mais próximo = NPNT (102,917 bp).

O5′ anotado de lncRNA LOC119872513 fica a32.785 bp de w01; o5′ de LOC102156669 a30.888 bp de w11. São extremidades anotadas, não TSS medidos experimentalmente. A ausência de overlap lncRNA não prova isolamento regulatório.

Na revisão [Ahmed2026, secção2](https://pmc.ncbi.nlm.nih.gov/articles/PMC12785581/), o critério referido de50 kb usa o5′ de genes codificantes. Portanto, a presença desses lncRNAs próximos **não constitui automaticamente falha desse critério**. A checklist local confirma que o código histórico extraiu genes codificantes, embora a palavra “ALL” fosse usada de forma ambígua. Deve separar-se esse critério de proximidade de RNA não codificante. Um [estudo primário de2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10836832/) aplicou, entre outros filtros, distância superior a100 kb de lncRNA; isso é uma regra mais restritiva distinta, não um limiar que se possa atribuir silenciosamente ao nosso pipeline. Sob essa regra de sensibilidade, ambas as janelas atuais falhariam. Não se mudou o filtro para preservar ou excluir candidatos nesta etapa.

## Trabalho pendente ativo

1. Integração das variantes PASS com unicidade exata concluída; ver atualização abaixo. AF por alelo normalizado concluída; ver atualização seguinte.
2. Avaliar especificidade com mismatches/PAM quando a nuclease estiver definida; foi perguntado ao utilizador qual nuclease e se existem sequências/genótipos dos cães.
3. Rever alcance/proveniência do contexto regulatório e SV, corrigindo afirmações excessivas, sem interpretar ausência de dados como aprovação.

O fecho experimental continua dependente do material individual e de evidência nas células relevantes. A automação de updates permanece **ativa**; não foi pausada pela produção deste checkpoint.

## Atualização: sensibilidade do contexto TAD e integração de máscaras

Recalculados intervalos entre fronteiras separadamente nas escalas 100, 250 e 500 kb, além da união histórica. Em w01 e w11, nenhum dos oito intervalos contém genes do catálogo de risco usado. Este resultado é robusto às três escalas testadas, mas continua a ser uma proxy: não demonstra isolamento tridimensional em células T caninas. O código de origem descreve os dados Hi-C como destinados a scaffolding de assembly; a identidade celular relevante não fica estabelecida por esta análise.

| Janela | Escala | Intervalo ROS (0-based, end-exclusive) | Genes de risco no catálogo |
|---|---|---|---|
| w01 | union | NC_051811.1:17125000-17350000 | nenhum |
| w01 | is_boundary_100000 | NC_051811.1:17050000-17350000 | nenhum |
| w01 | is_boundary_250000 | NC_051811.1:17125000-17650000 | nenhum |
| w01 | is_boundary_500000 | NC_051811.1:17025000-17675000 | nenhum |
| w11 | union | NC_051836.1:26900000-27350000 | nenhum |
| w11 | is_boundary_100000 | NC_051836.1:26900000-27350000 | nenhum |
| w11 | is_boundary_250000 | NC_051836.1:26900000-27450000 | nenhum |
| w11 | is_boundary_500000 | NC_051836.1:26800000-27475000 | nenhum |

A checklist histórica de 05_SHIP foi marcada como desatualizada. A conclusão “7/8 plenamente satisfeitos” não deve ser usada: o soft-score ATAC não valida acessibilidade em CAR-T, os intervalos entre fronteiras não validam TADs em CAR-T e já existem análises de conservação posteriores à checklist. Nenhum veto foi relaxado e as contagens anteriores foram preservadas.

Integração das duas máscaras já concluída em twenty_base_sequence_constraints.tsv: dos 981 segmentos de 20 bp por janela, 567 em w01 e 691 em w11 não têm ocorrência exata adicional nem sobreposição de REF de variantes PASS catalogadas. Isto não corresponde a guias selecionados, ausência de variantes no cão individual ou análise de mismatches/PAM.

AF específica concluída na atualização seguinte. Ainda pendente: revisão adicional do SV de w11; especificidade dependente da nuclease e dados individuais. A acessibilidade e segurança funcional em CAR-T continuam por demonstrar. Nenhum worker pesado está ativo neste checkpoint; esta etapa consistiu em auditoria local concluída.

## Frequências por alelo verificadas — 19 setembro, atualização seguinte

Concluída correspondência única dos 48 alelos normalizados aos ALT originais por equivalência do haplótipo local. A nova tabela `normalized_variant_alleles_with_AF_ROS.tsv` preserva o máximo por registo, acrescentando AF específica, índice do ALT original e identidade do alelo na referência UU. As tabelas anteriores permanecem intactas.

- w01: 27 alelos totais, 25 PASS, 10 alelos PASS com AF >=1%. O ALT T de w01_1 tem AF 0,302%, enquanto o máximo do registo é 13,5%; usar o máximo neste alelo teria alterado incorretamente a sua classe >=1%.
- w11: 21 alelos totais, 19 PASS, 9 alelos PASS com AF >=1%. O ALT TTGCTTC de w11_19 tem AF 0,05033%, versus máximo do registo 0,7801%; ambos abaixo de 1%.

Não confundir contagens de alelos com contagens de registos anteriores. Nenhuma variante foi retirada. As máscaras de REF de todas as variantes PASS e os resultados de unicidade exata permanecem iguais, porque não dependiam do limiar de AF. As AF são do catálogo, não dos cães que serão usados.

Revisitada a evidência já existente sobre o SV de w11: os quatro portadores catalogados têm 10–26 SNPs heterozigóticos de alta qualidade nos pequenos intervalos amostrados dentro da deleção alegada de 28,9 Mb. Isso questiona uma perda constitutiva simples nesses locais, mas a amostragem de SNPs não mede copy number de forma não enviesada, não verifica breakpoints e não exclui rearranjos complexos/mosaicismo. Não há nova validação por reads ou ensaio ortogonal nesta atualização. O SV mantém-se como limitação não resolvida; não foi convertido em falso positivo.

Próximas dependências: definir nuclease para especificidade com PAM/mismatches; localizar dados individuais se existirem; determinar se existe evidência bruta adicional adequada para o SV. A acessibilidade e função em CAR-T continuam a exigir dados relevantes. Sem worker pesado ativo nesta atualização; acompanhamento continua ativo.

## Reprodutibilidade e disponibilidade de dados — 19 setembro, 17:47 UTC checkpoint

Recomputadas independentemente as 1.962 linhas de variantes PASS versus unicidade exata de 20 bp: concordância integral com a tabela anterior. Script e recibo em verify_constraints_and_inputs.py e constraints_reproduction_and_input_inventory.json. Nenhuma nuclease escolhida foi encontrada na pesquisa de markdown/yaml/json do projeto; não se identificaram genótipos dos futuros dadores. A pesquisa de ficheiros no projeto não encontrou BAM/CRAM, sem afirmar ausência noutros discos.

Confirmados no ENA os quatro runs WGS paired dos portadores catalogados do SV: SRR15100208 (28,01 GB), SRR12330331 (32,60 GB), SRR15734832 (31,32 GB), SRR14750332 (36,17 GB), tamanhos decimais comprimidos. Total cerca de 128,1 GB: descarregar todos de uma vez excederia a folga atual em D. As respostas não fornecem links submitted_ftp de alinhamentos; isso não prova inexistência noutro repositório. Cada resposta individual ficou guardada; a consulta inicial com acessões separadas por vírgula retornou apenas cabeçalho e não foi tratada como ausência de dados. As consultas individuais confirmaram os runs e BioSamples esperados.

Não houve download de reads nem nova análise de copy number/breakpoints. Existe agora um caminho público concreto para investigar o SV, mas é uma etapa WGS substancial, que exige processamento de uma amostra de cada vez, espaço reservado e referência coerente. Estes portadores do catálogo não são os futuros dadores CAR-T. A próxima etapa útil é verificar ferramentas/referência/possibilidade de acesso regional antes de qualquer download grande. Acompanhamento continua ativo.

## Pré-verificação WGS/SV — 19 setembro, 17:55 UTC checkpoint

BWA, samtools 1.24 e os ficheiros de índice BWA de UU estão presentes. Ainda falta verificar integridade completa e correspondência da referência antes de novo alinhamento. Nenhum worker pesado ativo e nenhum download de reads iniciado.

As três listagens públicas de cram-share Dog10K foram inspecionadas; não incluem os quatro portadores em causa. Isto limita a via regional encontrada, sem provar inexistência de alinhamentos noutro local. Fonte: https://kiddlabshare.med.umich.edu/dog10K/cram-share/ .

| Portador | FT genótipo | GQ | DP | AD ref/alt |
|---|---|---|---|---|
| BFDB000001 | FAIL3 | 9 | 18 | [15, 3] |
| GBGV000001 | PASS | 42 | 7 | [4, 3] |
| PBGV000010 | PASS | 93 | 22 | [14, 8] |
| PBGV000012 | PASS | 87 | 19 | [6, 11] |

O FILTER PASS do registo não significa que todos os genótipos sejam PASS. BFDB000001 tem FAIL3/GQ9; portanto não deve ser escolhido primeiro só por ter o menor FASTQ. Prioridade técnica por FT e GQ: PBGV000010, PBGV000012, GBGV000001. Esta escolha prioriza suporte da chamada, não representatividade populacional. A chamada permanece não resolvida; a evidência de heterozigotia dispersa não foi promovida a refutação. Próximo passo: preparar análise limitada de um portador PASS, com plano de espaço e preservação de inputs não validados.

## Execução ativa: download de um portador PASS — 19 setembro, 18:04 UTC

Iniciado worker sv_download_one_carrier.py (PID inicial 6720) para PBGV000010 / SRR15734832 / SAMN21036425, selecionado por FT PASS e GQ93. São dois FASTQ comprimidos, total 31.315.288.945 bytes. Apenas esta amostra foi autorizada pelo worker; os mates são sequenciais. Estado em SV_PBGV000010/status.json e erros em SV_download_stderr.log. Usa lock exclusivo do sistema operativo, verificação de espaço com reserva mínima de 45 GiB, ficheiros .part preservados em falha e validação de tamanho/MD5 ENA antes de renomear. Nenhum input anterior é apagado. Não iniciar outro worker se o lock/processo estiver ativo.

Referência UU: SHA-256 795fbe3a51c4f48f0a08f475834a2836cad8c7690c70182aa1c7d594c1d7190a; 2.481.983.352 bases, 2.197 contigs. Todos os nomes/comprimentos FAI conferem com BWA .ann. Isto documenta consistência de metadados, não valida byte a byte o BWT nem substitui checksum original do fornecedor. Recibo UU_reference_preflight.json.

Ainda não há alinhamento ou nova evidência sobre o SV. Próxima etapa: acompanhar download e MD5; preparar análise do genoma completo com retenção de métricas/reads regionais, evitando o viés de alinhar apenas contra chr32. Não usar subset de reads como refutação de SV sem poder de cobertura demonstrado. Quando ambos os mates forem validados, dimensionar espaço de alinhamento antes de lançar. Acompanhamento de cinco minutos permanece ativo; w01/w11 continuam candidatos sem validação biológica.

## Monitorização — 19 setembro, 18:12 UTC

Worker único PID6720 confirmado ativo; sem erro no stderr. Primeiro mate tinha 511.705.088 bytes às 18:12:39 UTC, de 14.802.392.696 esperados. Tamanho confirmado abrindo o ficheiro e consultando o fim: o tamanho zero mostrado inicialmente pela listagem Windows era metadado temporariamente desatualizado do ficheiro aberto, não ausência de escrita. Zero mates validados até este checkpoint. Cerca de 119,38 GiB livres em D e 55,75 GiB em C. WSL disponibiliza cerca de 31 GiB RAM e 12 CPUs lógicos.

O débito observado é da ordem de 1 MB/s; os 31,3 GB totais podem demorar várias horas se este ritmo persistir. Não se promete conclusão iminente. O download permanece ativo e não foi duplicado.

Guardado SV_analysis_execution_design.json: alinhamento deverá usar o assembly inteiro, mesmo que só se retenham alinhamentos de chr32 e mates relevantes; cobertura distribuída e evidência de breakpoints devem ser avaliadas separadamente, com sensibilidade a duplicados/MAPQ e controlos de GC/mapeabilidade. Ainda não é um pipeline executado. O resultado de um portador do catálogo não substitui genótipo do dador CAR-T. Estado dos candidatos inalterado. Updates continuam ativos.

## Monitorização e preparação de QC — 19 setembro, 18:19 UTC

Transferência progrediu para 1.153.433.600 bytes no primeiro mate (7,79% dos 14.802.392.696 bytes; 3,68% dos 31.315.288.945 bytes totais). Um único worker PID6720, sem erros registados; 118,76 GiB livres em D. Nenhum mate validado e nenhum alinhamento iniciado. A estimativa continua na ordem de várias horas, dependente do débito.

Preparado validate_sv_fastq_pairs.py para verificar ambos os gzip até EOF, estrutura FASTQ, sincronização dos identificadores de pares e métricas descritivas de bases/N/Q30. Só aceita status download_complete e recibos de tamanho/MD5 do download; não executar nos .part. Testes sintéticos passaram para um par válido, rejeição de IDs desencontrados e rejeição de comprimento de qualidade incorreto. O script não foi executado nos dados reais, ainda incompletos. Essa validação é de integridade de input, não confirma identidade biológica nem resolve o SV. Manter o worker e acompanhamento ativos.

## Primeiro mate validado; espaço reavaliado — 19 setembro, 22:33 UTC

R1 completo: 14.802.392.696 bytes e MD5 7599c82bf5276df8e650cdeb9d05f730 conferem com ENA; recibo presente e tamanho do ficheiro final verificado. R2 em transferência, 4.731.174.912 de 16.512.896.249 bytes no checkpoint. Total recebido 19,53 GB, cerca de 62,4%. Worker único PID6720 ativo, stderr sem erros. Nenhum alinhamento nem validação FASTQ pareada ainda executados.

D: 63,29 GiB livres; C: 53,75 GiB. A redução em D desde o checkpoint anterior excede o aumento esperado destes FASTQ; a causa adicional não foi identificada e não foi atribuída ao worker. Restam cerca de 10,97 GiB para R2, projetando cerca de 52,3 GiB livres se não houver outro consumo. A reserva de 45 GiB do download continua ativa. Não lançar alinhamento antes de nova medição e dimensionamento dos temporários. Inputs não validados permanecem preservados.

Houve falha na comunicação dos updates entre os últimos heartbeats; não se pressupõe que tenham sido feitas verificações durante esse intervalo. Esta atualização baseia-se nos ficheiros/processo observados agora. Classificação w01/w11 inalterada, SV não resolvido. Próximo passo: concluir R2, verificar MD5 e validar os pares; acompanhamento permanece ativo.

## Preparação independente durante download — 19 setembro, 23:04 UTC

Download R2 ativo, 5,75 GB recebidos no início da etapa; R1 já validado. Preparadas 403 janelas consecutivas de até 100 kb cobrindo integralmente os 40.225.481 bp de NC_049253.1 (chr32 UU), sem gaps nem sobreposição. Destas, 288 ficam totalmente dentro do intervalo SV catalogado, 113 totalmente fora e duas cruzam os limites. Calculados GC entre bases canónicas, fração canónica e número de N por janela. Preparadas duas regiões de +/-10 kb em torno dos limites catalogados, preservando as coordenadas 0-based originais.

Ficheiros: SV_chr32_100kb_reference_windows.tsv, SV_chr32_100kb_windows.bed, SV_breakpoint_neighborhoods_UU.bed e SV_coverage_window_preparation.json. Script prepare_sv_coverage_windows.py. Esta é preparação de referência, não cobertura medida da amostra. Janelas externas não são controlos diploides validados; falta avaliar mapeabilidade/repeats e a cobertura real. Nenhuma conclusão SV ou de safe harbor foi alterada. A análise continua dependente da conclusão e validação de R2; espaço será reavaliado antes de alinhar.

## Recuperação de ligação — 19 setembro, 23:11 UTC

Worker anterior terminou às 23:06:13 UTC por ConnectionResetError WinError10054 (ligação fechada remotamente). Não houve erro de disco ou MD5. R1 validado preservado; R2 parcial preservado, último progresso registado 5.838.471.168 bytes. Recibo de falha arquivado em SV_PBGV000010/failure_20260919_230613.json; logs originais mantidos.

Confirmada ausência do worker antes de relançar. Novo PID16284 usa o mesmo lock exclusivo e código de retoma por HTTP Range, exigindo resposta206 e Content-Range coerente antes de acrescentar bytes. Revalida o R1 já existente antes de retomar R2. Logs novos: SV_download_resume1_stdout.log / SV_download_resume1_stderr.log. D com62,29GiB livres; reserva45GiB permanece. Não se declara retoma confirmada até observar progresso do novo PID no status.json. Nenhuma alteração biológica da classificação; falta concluir download e análise do SV.

## Segunda recuperação da ligação — 20 setembro, 07:08 UTC (heartbeat recebido com timestamp00:25)

Nova interrupção remota WinError10054 às00:21:25UTC, com último progresso de R2 em13.300.137.984bytes. Confirmada ausência de worker anterior antes de relançar PID7856 com lock exclusivo. Falha preservada em SV_PBGV000010/failure_20260920_002125.json; logs anteriores intactos; logs novos SV_download_resume2_stdout.log e SV_download_resume2_stderr.log. R1 permanece preservado e é revalidado antes da retoma por Range. Nenhum dado parcial apagado. D:55,34GiB livres. A queda invalida a previsão anterior de tempo; não há compromisso de hora de conclusão. Próximo checkpoint deve confirmar status do novoPID e avanço deR2, depois tamanho/MD5 e validação pareada. Nenhum alinhamento ou nova conclusão biológica nesta etapa.

Correção temporal: relógio real confirmado em07:08UTC/09:08local. A retoma foi efetuada agora, não às00:25. Status novo PID7856 confirma download retomado às07:08:23UTC, R2=13.304.332.288bytes. O período sem execução não deve ser relatado como monitorização contínua.

## Retoma após timeout — 20 setembro, 07:54 UTC

Worker anterior parou às07:48:13UTC por timeout de leitura, R2 último progresso15.955.132.416bytes. Processo ausente confirmado. Relançado PID27164 com lock exclusivo e retoma Range; ficheiro parcial e falha failure_20260920_074813.json preservados. Logs atuais SV_download_resume3_stdout.log / SV_download_resume3_stderr.log. R1 é revalidado antes da retoma. D52,87GiB livres. Ainda não há MD5 final de R2 nem validação pareada. Não declarar download completo até recibos. SV e classificação dos candidatos inalterados.

## Download concluído; validação pareada ativa — 20 setembro, 08:13 UTC

Ambos os mates concluídos e validados por tamanho e MD5 ENA às08:12:59UTC: R1 14.802.392.696bytes/7599c82bf5276df8e650cdeb9d05f730; R2 16.512.896.249bytes/597f33fb4418c2c741f2be55cc99baea. Total31.315.288.945bytes. Nenhuma cópia descomprimida criada, nenhum input apagado.

Iniciada validação integral FASTQ pareada PID26168, script validate_sv_fastq_pairs.py, lock exclusivo pair_validation.lock. Estado em SV_PBGV000010/pair_validation_status.json; logs SV_pair_validation_stdout.log e SV_pair_validation_stderr.log. Progresso cada100.000pares; verifica estrutura, IDs sincronizados, gzip até EOF e métricas de bases/N/Q30. Otimização da contagem Q30 verificada em teste sintético, incluindo rejeição de IDs desencontrados. Não confundir processo iniciado com validação concluída; receipt final fastq_pair_validation.json ainda pendente.

D tinha cerca52,3GiB livres. Alinhamento ainda não iniciado e exige reavaliação do espaço/temporários. Os403intervalos de referência e neighborhoods de breakpoints estão preparados. Nenhum novo resultado biológico sobre w01/w11 ou SV. Acompanhamento permanece ativo.

## Via prioritária w01 — 20 setembro

Preparado W01_TEST_HANDOFF/LEIA_PRIMEIRO.md, com27alelos normalizados e567segmentos de20bp reunidos em22grupos de inícios consecutivos. São segmentos sem flags nas duas verificações limitadas, não guias selecionados. Dossier distingue locus de ponto de corte, sangue bulk de célulasT, e define dependências concretas de nuclease, cassete e sequência dos futuros dadores. w01 avança em preparação independentemente doSVde w11. Validação FASTQde w11 continua ativa, sem novo resultadoSV.

## FASTQ concluído e alinhamento iniciado — 20 setembro, 12:14 UTC

Validação integral terminada às09:01:52UTC: 185,353,881 pares, 27,803,082,150 bases por mate, gzip atéEOF e identidade dos pares verificados. Q30:R1=94.50%;R2=87.76%. Isto é QC técnico, não confirmação de identidade biológica. Houve intervalo sem updates/lançamento entre a conclusão do QC e esta execução; não relatar como alinhamento contínuo.

Lançado agora align_sv_one_carrier.py, supervisor WSL PID294, launcher Windows26044, lock exclusivo. Alinhamento BWA do par completo contra todo o assembly UU; retenção de reads emchr32 ou com mate emchr32. Estado SV_PBGV000010/alignment/status.json; logs pipeline_0/1/2.stderr.log e SV_alignment_launcher_stderr.log. Teste sintético confirmou retenção de ambos os mates discordantes e exclusão do par não relacionado, além de sort/fixmate/markdup/index. BAM final ainda não existe.

Reserva: exige30GiB livres no arranque e interrompe abaixo de20GiB durante execução, preservando todos os inputs/intermediários. Isto é política própria do alinhamento regional; o download usava45GiB. Livre no arranque cerca52,35GiB. Nome-sort, fixmate, coordinate-sort, markdup, quickcheck/index seguem automaticamente apenas após sucesso da etapa anterior. Duplicados são marcados, não removidos. A marcação no subset retido não equivale à marcação do genoma completo para todos os fragmentos inter-cromossómicos; manter essa limitação na revisão. Nenhuma evidência de cobertura/breakpoint interpretada ainda. w01 continua independente desta investigação de w11.

## Análise pós-alinhamento preparada — 20 setembro

Script summarize_sv_alignment.py preparado, ainda não executado em dados reais. Requer completed.json e BAM indexado. Calcula cobertura de bases alinhadas nas403janelas para MAPQ20/30, incluindo/excluindo duplicados marcados; mantém janelas parciais e N explícitos, e apresenta separadamente resumo descritivo de janelas completas com>=99%bases canónicas. Contagens de flags SA/supplementary/paired-not-proper nas vizinhanças de breakpoints servem apenas para revisão. Sem correção GC/mapeabilidade e sem inferência automática de deleção; mates sobrepostos contam duas vezes. Testes da distribuição de blocos entre janelas e coordenadas half-open passaram. Nenhum novo resultado biológico. Executar apenas após conclusão do alinhamento e rever antes de interpretar.

## Falha SAM e diagnóstico ativo — 2026-09-20T17:21:11.691682+00:00

Alinhamento original falhou às15:55:10UTC, samtools view reportou Parse error at line219236384. BWA tinha registado218.133.606reads processados (58,84%); isto não representa alinhamento final validado. Sem processos antigos ativos ao verificar. Nenhum indício de falta de disco nos logs; causa exata permanece desconhecida. O BAM parcial e temporários não foram apagados nem aceites como resultado.

Iniciado diagnose_sv_sam_failure.py: extrai pares106.000.000–109.999.999 (índices0based) de ambos os FASTQ validados, repete BWA contra todoUU e guarda SAM diagnóstico, depois testa parsing comsamtools. O intervalo é uma aproximação ao lote que falhou, não identificação exata do read problemático; diferenças de batch podem afetar reprodução. Se o erro reproduzir, guarda contexto da linha. Piso20GiB de espaço durante execução, lock exclusivo. Estado SV_PBGV000010/sam_failure_diagnostic/status.json; logs SV_sam_diagnostic_stderr.log. Não relançado o alinhamento completo nem duplicado processo.

Falharam novamente os updates regulares; não se declara acompanhamento contínuo no intervalo. w01 continua preparação independente, nuclease/cassete/genótipo pendentes. w11SV continua não resolvido.

## Diagnóstico SAM — 2026-09-20 18:35 UTC
Repetição limitada concluída: 4.000.000 pares (índices zero-based 106.000.000–109.999.999), alinhados ao genoma UU completo. Conversão SAM→BAM terminou com exit 0; samtools quickcheck exit 0; leitura integral por samtools view -c -F 2304 contou exatamente 8.000.000 registos primários. BAM de diagnóstico: 863.267.650 bytes. Não reproduziu a falha nesta configuração. Não demonstra a causa original nem valida os outputs parciais da execução completa. Ambos os processos de diagnóstico terminaram; nenhum alinhamento completo de recuperação foi iniciado. Próximo passo: recuperação com blocos e recibos verificáveis, preservando inputs e outputs anteriores. w11 SV permanece não resolvido; w01 independente. Espaço observado: D 47,15 GiB; C 43,36 GiB.

## Recuperação por blocos iniciada — 2026-09-20
Worker recover_sv_chunks.py ativo em SV_PBGV000010/alignment_chunked. 185.353.881 pares em 47 blocos, até 4 milhões por bloco. Cada bloco alinha contra todo o genoma UU, guarda SAM temporário, valida leitura integral e contagens primárias/R1/R2, retém alinhamentos chr32 ou mate chr32 e grava BAM com SHA256/recibo. Só os FASTQ e SAM gerados de blocos validados são removidos; originais, diagnóstico e tentativa falhada preservados. Todos os blocos serão reunidos antes de fixmate e markdup. Estimação de insert-size por bloco pode diferir da corrida monolítica. Teste sintético end-to-end passou: 2 blocos, 4 registos primários finais, indexação e preservação dos originais. Retoma valida hashes e percorre fontes até ao último bloco validado; um bloco parcial não validado exige investigação e não é sobrescrito. Ainda não existe alinhamento completo validado nem conclusão SV. Não duplicar; consultar status.json antes de qualquer lançamento.


## Retoma após reinício — 2026-09-21T07:45:15.3734225Z
Confirmados 19 recibos e hashes dos BAMs, 76.000.000 pares validados. Nenhum worker anterior ativo. chunk_019 tinha apenas dois FASTQ gerados parciais (~78 MB cada), preservados em interrupted_chunk_019 com timestamp. Worker relançado com lock exclusivo. Percorre os FASTQ originais e salta os blocos validados antes de reextrair chunk_019; não realinha os 19 concluídos. Espaço no arranque D46,75 GiB/C47,66 GiB. SV ainda sem conclusão. Monitorização mantém frequência horária.



## Retoma após shutdown — 2026-09-21T12:07:07.3244506Z
23 blocos/92.000.000 pares: hashes e contagens conferidos. Nenhum worker anterior ativo. chunk_023 interrompido durante alinhamento, preservado integralmente em interrupted_chunk_023 com timestamp (dois FASTQ de 1.348.000.000 bytes e SAM parcial de 2.997.735.424 bytes). Worker relançado; percorre originais até ao ponto de retoma, sem realinhar blocos validados. Espaço D41,40 GiB/C52,71 GiB. Sem conclusão biológica nova; acompanhamento horário mantido.


## Revisão GC e mapeabilidade — 2026-09-22 15:34 UTC
39.997 probes canónicas de 150bp remapeadas ao assembly completo. Controlos GC +/-0,01 e proxy de mapeabilidade +/-0,05: 281/288 janelas internas com 3–5 externas; razão mediana 0,99806. Verificação independente das correspondências, razões e hashes passou. 7 janelas sem controlos preservadas como não emparelhadas. Isto conclui esta comparação descritiva WGS, não o painel experimental de controlos nem a auditoria global. Relatório W11_WGS_REVIEW_2026-09-22.md reúne evidência, limitações e próximos passos. w11 condicional; w01 prioritário. Sem worker pesado ativo.

