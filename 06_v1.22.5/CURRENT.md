# Safe Harbor CAR-T — v1.22.5: desenho emparelhado, genes ausentes e dados brutos

16 de setembro de 2026. **O artigo confirma três cães com condições emparelhadas. A pesquisa oficial de sinónimos não recuperou os 15 genes ausentes da matriz processada. Os nove runs de RNA-seq bruto foram localizados e associados às amostras.** Esta etapa esclarece a origem dos dados; não modifica os scores ou valida candidatos.

## Desenho confirmado

O texto e a legenda da figura 6 do [artigo original](https://doi.org/10.1007/s00262-024-03642-4) descrevem RNA-seq de três cães e comparação emparelhada de CAR-T B7-H3, CAR-T B7-H3/CXCR2 e T não transduzidas. Assim, o desenho é três indivíduos × três condições, não nove indivíduos. O método específico de RNA indica colheita após dez dias de expansão após transdução; a descrição geral das culturas admite dias 10–14. Esta distinção fica agora registada.

A correspondência de B/E/M aos rótulos de réplica do GEO foi verificada na v1.22.4. A confirmação do desenho no artigo sustenta a interpretação emparelhada ao nível do estudo. Não prova a identidade física de cada tubo. Os quatro campos source/tissue contraditórios continuam explícitos; não foram corrigidos por suposição ou escondidos.

## Auditoria de nomes génicos

Foi obtida a tabela corrente [NCBI Gene para Canis familiaris](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Canis_familiaris.gene_info.gz). Para os 40 genes de contexto foram pesquisados símbolo, símbolo oficial de nomenclatura, sinónimos e identificador LOC associado ao GeneID. Só seria aceite recuperação automática com correspondência única de gene e de linha da matriz.

**Nenhum dos 15 genes ausentes foi recuperado por um sinónimo inequívoco na matriz CPM.** Entre eles permanecem ANO2, NTF3, CD247 e CD8B. O resultado fecha a hipótese simples de recuperar esses valores apenas renomeando linhas com os sinónimos correntes consultados; não determina se a ausência veio de anotação histórica, filtragem ou baixa expressão. A tabela NCBI atual não é necessariamente a versão de anotação dos autores.

As matrizes anteriores mantêm esses valores ausentes, não zero. Não foi inferido silêncio génico nem ausência de função regulatória.

## Disponibilidade real de dados

O diretório público de suplementos de GSE247355 consultado contém apenas a matriz CPM já usada. O material suplementar ligado ao artigo é um PDF de figuras; não identificámos ali uma matriz de contagens brutas. Os métodos referem STAR/CanFam3.1, HTSeq e processamento RUVseq/DESeq2, mas não fornecem nessa descrição o ficheiro/release exato de anotação ou a regra que produziu exatamente as 14.385 linhas da matriz disponibilizada. Não se pode prometer reprodução exata a partir destes metadados.

Os nove SRX foram resolvidos no [ENA](https://www.ebi.ac.uk/ena/browser/home), cada um com um SRR e dois FASTQs. A associação GSM→SRX→SRR/BioSample foi verificada, com sample_alias correspondente ao GSM. Os ficheiros somam **63.892.132.810 bytes, aproximadamente 63,9 GB decimais comprimidos**. O manifesto tem condições, rótulos emparelhados, accessions, URLs e tamanhos. Estes aliases confirmam ligações de registos, não resolvem independentemente os campos biológicos contraditórios.

Os FASTQs não foram descarregados nem reprocessados nesta etapa. Reconstruir contagens a partir deles é uma nova análise computacional substancial, não uma conversão algébrica válida de CPM para contagens inteiras. A matriz CPM não foi utilizada em testes que exigem contagens brutas.

## Trabalho concretamente preparado

`raw_RNA_reprocessing_manifest.tsv` é a entrada rastreável para eventual reprocessamento dos nove runs. Será necessário fixar referência/anotação e versões, verificar integridade e QC, produzir contagens por gene e rever identidade amostral antes de usar um modelo emparelhado. Uma análise com anotação ROS atual deverá ser rotulada como reanálise, não reprodução idêntica do pipeline CanFam3.1 original.

A utilidade seria recuperar contexto de genes ausentes e comparar condições com contagens apropriadas. Continuaria a não medir diretamente acessibilidade das janelas nem neutralidade de integração. A revisão atual dos candidatos permanece exploratória e a necessidade de validação em T/CAR-T relevantes não foi removida.

## Ficheiros e verificação

`gene_alias_audit.tsv`, `source_audit_summary.json`, `ENA_raw_read_inventory.json`, nove respostas ENA TSV, `raw_RNA_reprocessing_manifest.tsv`, `paper_method_evidence.json`, texto XML do artigo, inventários dos diretórios públicos e snapshot NCBI Gene com hashes. As nove associações e pares FASTQ foram verificados; nenhuma recuperação de gene foi inventada. Scoring global permanece v1.21.8; contexto descritivo RNA permanece v1.22.4 com desenho agora confirmado no artigo.
