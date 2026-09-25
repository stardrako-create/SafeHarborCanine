# w01 — preparação para teste, não validação de safe harbor

**Decisão:** preparar w01 como primeira prioridade experimental, em paralelo à investigação separada do SV de w11. O SV de w11 não constitui evidência contra w01 e a sua análise não é pré-requisito para preparar w01. Nenhum critério foi relaxado.

## Locus e identidade

- ROS_Cfam_1.0: NC_051811.1:17255706–17256706, BED 0-based/end-exclusive.
- A mesma janela em coordenadas 1-based/inclusivas: NC_051811.1:17255707–17256706.
- UU_Cfam_GSD_1.0: NC_049228.1:17403010–17404010, BED, orientação positiva; correspondência integral de sequência de 1 kb verificada.
- A janela é a região a investigar, não um ponto de corte escolhido.

## Evidência já disponível

ATAC de sangue: percentil 63,75 comparado com janelas de tamanho equivalente. Não equivale a percentil em células T. Dois inícios de fragmentos em dois núcleos T no conjunto exploratório de um único dador; isto não demonstra enriquecimento robusto. RRBS: apenas20% da janela observada, metilação média observada60,71%; não designar como baixa metilação comprovada. Sem overlap de repeats no teste atual. Contexto de transcritos inclui lncRNA próximo: ausência de overlap não é isolamento regulatório.

Variantes:27 alelos normalizados,25 PASS,10 alelos PASS comAF>=1%. As frequências são populacionais do catálogo. O mapa anexo reúne567 segmentos consecutivos de20bp com uma ocorrência exata no assembly e sem sobreposição de REF de variantesPASS catalogadas. Não são567 guias: falta nuclease/PAM, mismatches/bulges, contexto dos indels e genótipo individual. Intervalos de starts não autorizam qualquer corte em toda a sua união.

## O que falta para um desenho de edição executável

1. Definir nuclease/variante e estratégia de integração; não foi encontrada decisão explícita no projeto.
2. Definir o cargo/cassete e objetivo de expressão. Sem isso não se fecha desenho do donor nem ensaios específicos da integração.
3. Confirmar a sequência e variantes do locus no material dos cães a usar, incluindo contexto necessário ao desenho; registar referência e convenção de coordenadas.
4. Escolher e verificar o alvo dependente da nuclease: PAM, semelhanças genómicas, variantes nos elementos do desenho e especificidade das sequências dos ensaios.

## O que a experiência precisa de estabelecer

Comparar células editadas com controlos adequados para verificar identidade e estrutura da integração, persistência e heterogeneidade da expressão, efeito nos genes próximos e manutenção de viabilidade e função das células T/CAR-T. Os critérios de aceitação devem ser definidos antes de examinar os resultados, em função da cassete e aplicação. Acessibilidade/atividade em células T relevantes é uma lacuna real, não preenchida pelo sangue bulk. Este documento não contém parâmetros operacionais de edição, guias ou primers validados.

## Fontes internas

07_FINAL_CANDIDATES_2026-09-19/shortlist_evidence.tsv;08_CANDIDATE_CLOSURE_AUDIT/variant_regulatory_review.json;allele_frequency_review.json;exact_sequence_review.json;coding_vs_all_gene_context.json;tad_scale_audit.json. O relatório histórico do painel não representa fecho biológico.
