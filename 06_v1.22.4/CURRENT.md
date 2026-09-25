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
