# v1.21.2 — reconciling two parallel working sessions

2026-09-12. Este projeto teve duas sessões de trabalho independentes a
correr sobre o mesmo código no mesmo dia: esta (Claude) e uma paralela
(Codex, ver `06_v1.21.1/`). Este ficheiro documenta o que cada uma
encontrou, para que nenhuma correção real se perca e nenhuma fique por
atribuir.

## O que esta sessão (Claude) fez, antes de saber da outra

- Corrigiu `veto_gene_dense_neighborhood` para medir distância à
  extremidade 5' dos genes (verificado diretamente contra o texto de
  Ahmed et al. 2026, Secção 2), não ao corpo do gene - `scripts/
  extract_genes_stranded.py` (novo) + `distance_to_nearest_5prime_end()`.
  380→264 exclusões por este veto.
- Construiu uma anotação de ilhas CpG própria e auditável
  (`scripts/build_cpg_islands.py`) para verificar de forma independente o
  BED regulatório opaco do Ehsan (75.600 intervalos sem tipo/evidência) -
  55% de sobreposição com o conjunto dele, validação real mas parcial (só
  cobre promotores, não enhancers distais).
- Redesenhou a análise de sensibilidade depois de um erro circular ser
  apanhado (pesos não podem alterar vetos duros; "100% robusto" com p55
  era um artefacto de campo vazio, não robustez real).

## O que a outra sessão (Codex) fez, documentado em 06_v1.21.1/

- **Achado real que esta sessão não tinha apanhado**: `rrbs_mean` incluía
  bases sem cobertura como zero medido (a mesma classe de bug já corrigida
  aqui para `atac_variability`/`peak_frequency`, mas nunca aplicada a
  `rrbs_mean` em si). Corrigido via `evidence_summary()` em `bw_utils.py`.
- Adicionou estados de evidência explícitos (`evaluation_status`/
  `missing_evidence`) - falta de anotação (repeat/conservation/TAD/etc.)
  deixa de contar silenciosamente como aprovação.
- Escreveu `06_v1.21.1/CAR_T_ANALYSIS_2026-09-12.md`, uma análise cuidada
  do objetivo CAR-T especificamente - levanta que os 71 cães são PBMC, não
  células T, entre outros pontos ainda por fechar. Vale a pena ler na
  íntegra antes de qualquer decisão biológica.
- Escreveu diretamente em MemPalace (wing "Locus Canino", rooms "Overview"
  e "CAR_T_Protocol_2026-09-12") - conteúdo já verificado, coerente com o
  código encontrado.

## Como foi reconciliado

O código de `06_v1.21.1/code/` já partia do código desta sessão (os campos
`gene_5prime_clearance` já lá estavam) - não são duas implementações
independentes do zero, é uma sessão a construir sobre a outra. As correções
próprias de Codex (`evidence_summary`, o guard de `attempts` em
`build_matched_window_background`, o `ValueError` explícito em intervalos
inválidos, `evaluation_status`/`missing_evidence`) foram revistas e
adotadas nos scripts principais (`scripts/bw_utils.py`,
`scripts/score_ship_candidates_v2.py`, `scripts/build_methylation_track.py`
- 5 testes de regressão novos, 11/11 a passar).

**Verificação cruzada**: correr o scoring principal já fundido contra as
mesmas tracks reproduziu exatamente os números que Codex reportou -
ANO2/NTF3 score 0.5072, mediana do background RRBS 77.15, os mesmos 3
candidatos em `insufficient_evidence` com os mesmos motivos. Confiança
real de que as duas linhas de trabalho convergem no mesmo resultado, não
apenas por coincidência de wording.

## O que fica por fazer, de ambas as sessões

- Repeat_content/conservation para os 3 candidatos novos (bloqueador para
  os tirar de `insufficient_evidence`).
- A auditoria de PBMC-vs-T-cells e geometria do ponto de inserção local
  proposta em `CAR_T_ANALYSIS_2026-09-12.md`, secções 2 e 6.
- Reconstruir as tracks ATAC com a recuperação de blocos corrompidos
  (código já corrigido nesta sessão, 2026-09-11 - tracks atuais ainda
  usam o fallback mais grosseiro, impacto já quantificado como negligenciável
  para RIT2 especificamente, não medido para o resto do genoma).
- Reconciliar as duas reconstruções de ilhas CpG independentes (113.459
  regiões desta sessão vs. 153.908 chamadas de Codex - convenções de
  contagem diferentes, não necessariamente um conflito, mas não verificado).
- Enviar o email ao Ehsan (rascunho existe, `05_SHIP/ehsan_reply_draft_
  2026-09-11.md` - desatualizado outra vez, precisa de refletir v1.21.2
  antes de enviar) - continua bloqueado no Chrome.
