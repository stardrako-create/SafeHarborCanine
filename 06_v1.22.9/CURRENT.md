# v1.22.9 — contexto CAR-T dos 40 genes previamente selecionados

19 setembro 2026. **Os 40 genes têm correspondência exata e única na anotação corrigida; os 15 ausentes da matriz publicada podem agora ser avaliados nesta reanálise.** Não são 15 genes recuperados pelo bug dos 60 genes: são problemas diferentes. A primeira ausência refere-se à tabela publicada; a segunda era a ingestão do GTF pela v1.22.6.

## Âmbito e validação

Mantida a seleção definida na v1.22.4: dez genes distintos mais próximos de cada uma das três âncoras declaradas, mais marcadores/genes predefinidos. Não se selecionaram genes pela expressão observada. Mantidas as âncoras históricas, sem as apresentar como coordenadas finais de integração. Os nomes históricos dos loci não substituem a lista efetiva de vizinhos e distâncias em 06_v1.22.4/local_gene_context.tsv.

Usada a matriz não orientada corrigida v1.22.7, cuja consistência foi verificada. Correspondência pelo gene_name único do índice ROS, sem aliases forçados. Produzidos 360 registos gene/amostra, 120 resumos gene/condição e 360 comparações descritivas dentro do mesmo rótulo de cão. Os totais foram verificados e os hashes de inputs guardados. A v1.22.8 já demonstrou melhor correspondência 9/9 com os rótulos publicados; isto não confirma identidade física.

CPM é contagem bruta × 1 milhão / soma dos pares atribuídos a genes na mesma biblioteca. É normalização descritiva por tamanho, sem correção de composição; não são valores TPM, nem necessariamente o CPM do autor. Comparações usam três cães reportados, B/E/M, cada um com controlo, CAR simples e CAR duplo.

## Resultados

| Gene | Ausente na tabela publicada | Amostras com contagem >0 | Contagens brutas min–max | CPM min–max |
|---|---|---:|---:|---:|
| AIMP1 | não | 9/9 | 1611–3342 | 41.9817–67.7073 |
| ANO2 | sim | 4/9 | 0–4 | 0.0000–0.0702 |
| APOBEC4 | sim | 6/9 | 0–6 | 0.0000–0.1072 |
| ARHGEF38 | sim | 8/9 | 0–6 | 0.0000–0.1237 |
| ARPC5 | não | 9/9 | 18112–30295 | 444.2237–579.3297 |
| C7H1orf21 | não | 9/9 | 755–1823 | 16.7493–33.1226 |
| CD247 | sim | 9/9 | 6800–16027 | 150.8550–322.9327 |
| CD3D | não | 9/9 | 4059–9139 | 90.0471–202.8930 |
| CD3E | não | 9/9 | 17770–40147 | 394.2197–830.3096 |
| CD3G | não | 9/9 | 4246–12170 | 94.1957–239.5205 |
| CD4 | não | 9/9 | 3609–13406 | 63.3619–297.6238 |
| CD8A | não | 9/9 | 2291–21131 | 50.8248–377.5567 |
| CD8B | sim | 9/9 | 469–17105 | 10.4046–305.6224 |
| CD9 | não | 9/9 | 19–247 | 0.3336–5.4836 |
| COLGALT2 | sim | 7/9 | 0–13 | 0.0000–0.2886 |
| EDEM3 | não | 9/9 | 2039–5525 | 54.8190–116.1322 |
| GALNT8 | sim | 5/9 | 0–5 | 0.0000–0.1075 |
| GIMD1 | sim | 9/9 | 1–8 | 0.0222–0.1776 |
| GSTCD | não | 9/9 | 152–979 | 4.0866–21.7187 |
| INTS12 | não | 9/9 | 476–1094 | 12.1250–19.2070 |
| KCNA1 | sim | 2/9 | 0–5 | 0.0000–0.0893 |
| KCNA5 | sim | 0/9 | 0–0 | 0.0000–0.0000 |
| KCNA6 | sim | 0/9 | 0–0 | 0.0000–0.0000 |
| LCK | não | 9/9 | 5653–28438 | 125.4093–508.1140 |
| LOC119867197 | sim | 9/9 | 18–53 | 0.4662–1.1140 |
| LOC119872513 | sim | 0/9 | 0–0 | 0.0000–0.0000 |
| LOC119876429 | sim | 8/9 | 0–32 | 0.0000–0.6990 |
| NCF2 | não | 9/9 | 28–1138 | 0.7083–24.8467 |
| NIBAN1 | não | 9/9 | 13471–22915 | 277.3597–481.6595 |
| NPNT | não | 7/9 | 0–13 | 0.0000–0.2886 |
| NTF3 | sim | 1/9 | 0–4 | 0.0000–0.0808 |
| PLEKHG6 | não | 9/9 | 8–75 | 0.1775–1.6400 |
| PPA2 | não | 9/9 | 1850–3340 | 38.8859–60.4111 |
| RGL1 | não | 9/9 | 7–279 | 0.1882–5.8644 |
| SMG7 | não | 9/9 | 6253–12915 | 168.1133–226.7439 |
| TBCK | não | 9/9 | 741–1187 | 18.1774–23.6933 |
| TET2 | não | 9/9 | 2371–3494 | 55.7942–67.0142 |
| TNFRSF1A | não | 9/9 | 276–1338 | 6.1229–25.0570 |
| TSEN15 | não | 9/9 | 1642–2746 | 37.1203–48.2105 |
| VWF | não | 5/9 | 0–8 | 0.0000–0.1650 |


CD247 e CD8B apresentam contagens em todas as amostras apesar de ausentes da tabela CPM publicada. Isso impede interpretar aquela ausência como ausência de expressão. ANO2 (0–4 pares) e NTF3 (0–4) têm suporte de contagem muito baixo; ratios relativos neste domínio são frágeis. LOC119876429 apresenta 0–32 pares, e LOC119872513 zero em todas as amostras. KCNA5 e KCNA6 também têm zero pares atribuídos. Contagem zero é uma observação sob esta biblioteca, anotação e regra de atribuição — não prova de silêncio absoluto.

## Comparações pareadas

within_donor_descriptive_ratios.tsv contém CAR simples/controlo, CAR duplo/controlo e duplo/simples por cão. O log2-ratio só existe quando ambos os CPM são estritamente positivos. Não se adicionou pseudocount para criar razões artificiais nas contagens zero. Valores ausentes por zero e por anotação têm estados explícitos. Os ratios de genes com pouquíssimos pares devem ser lidos com os respetivos numeradores/denominadores, não como efeitos fiáveis.

condition_descriptive_summary.tsv apresenta mediana e amplitude entre três cães em cada condição. Não foram calculados p-values, declarada estabilidade/equivalência, nem usados estes resumos como exclusão automática. Ausência de significância num eventual teste futuro também não provaria neutralidade.

## Implicações para os candidatos

Esta análise fecha a lacuna de contexto RNA dos 15 símbolos ausentes da tabela publicada. Não fecha o problema ATAC nem demonstra expressão estável de um transgene. Genes vizinhos com contagens substanciais continuam a exigir avaliação de possível perturbação regulatória; baixa expressão de outro vizinho não torna o intervalo automaticamente seguro.

Não foi alterado o ranking, nenhum filtro de candidatos foi afrouxado e nenhuma âncora foi promovida a alvo final. Continuam relevantes os conflitos regulatórios, cobertura e identidade celular do ATAC, incertezas de variantes populacionais e a falta de validação funcional pós-integração documentados nas auditorias anteriores.

## Pendências reais

Reconciliação do kit/protocolo e dos quatro conflitos source/tissue; correspondência génica mais ampla se for necessária inferência transcriptómica; normalização apropriada e modelo pareado se houver uma pergunta de expressão diferencial; evidência celular e funcional de acessibilidade, expressão persistente e segurança da integração. Nenhum contacto externo foi enviado nem nova execução pesada iniciada nesta etapa.

Artefactos: gene_mapping.tsv; all40_gene_sample_context.tsv; all40_overview.tsv; condition_descriptive_summary.tsv; within_donor_descriptive_ratios.tsv; context_summary.json. Scripts context_v1229.py e finish_context_v1229.py. A recontagem concluída continua a ser a v1.22.7.
