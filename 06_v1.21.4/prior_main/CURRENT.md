# Versão corrigida: v1.21.3

Execução concluída em 2026-09-12. v1.21.3 preenche repeat_content e
conservation (phyloP) para os 3 candidatos que v1.21.2 tinha deixado como
`insufficient_evidence` por falta dessas anotações - mesmo código de
scoring de v1.21.2, sem mudança de método. v1.21.2 juntou duas linhas de
trabalho paralelas sobre o mesmo projeto (Claude, nesta sessão; Codex,
sessão independente em `06_v1.21.1/`) - ver `06_v1.21.2/RECONCILIATION.md`
para a proveniência completa de cada correção anterior.

461 candidatos: 457 excluídos, **1 com evidência insuficiente** (loc2
`NC_051820.1:22976673-23048542`, LOC111090199/LOC100682550 - falta só
conservation: o intervalo não tem bloco ortólogo coerente em CanFam3,
fragmenta-se em 7 pedaços espalhados por 6 cromossomas ao liftOver; ver
`06_v1.21.3/CURRENT.md` para o detalhe - é discordância estrutural real,
não ausência de dados), **3 passam as verificações registadas**:
ANO2/NTF3 (`NC_051831.1:39677324-39739751`, score 0.5072, #1), NPNT/TBCK
(`NC_051836.1:27013911-27071204`, score 0.5015, #2), LOC119876429/
LOC119872513 (`NC_051811.1:17217706-17274511`, score 0.3988, #3). Os
scores não mudaram face a v1.21.2 - repeat_content/conservation são vetos
passa/falha, não entram na fórmula do score. **Nenhum é um safe harbor
validado.**

Código e resultados em `06_v1.21.3/`. Resultados legados (`06_v1.20.0/`,
`06_v1.21.0/`, `06_v1.21.1/`, `06_v1.21.2/`) mantidos intactos como
histórico - nunca sobrescritos.

Correções acumuladas até v1.21.2 (linhagem completa): janela vs. bin no
piso de acessibilidade; viés de amostragem de cromossomas no background;
distância aos 5' dos genes (Ahmed et al. 2026, não corpo do gene);
evidência insuficiente ≠ estabilidade (variabilidade ATAC/RRBS); média
RRBS não deve incluir bases sem cobertura como zero medido; estados de
evidência em falta agora explícitos (`evaluation_status`/
`missing_evidence`), não silenciosamente aprovados.

Por fazer antes de qualquer conclusão biológica: reconciliar a auditoria
dos 461→3 com a análise CAR-T (`06_v1.21.1/CAR_T_ANALYSIS_2026-09-12.md`)
que levanta pontos ainda por fechar (PBMC vs. T cells, geometria local do
ponto de inserção, auditoria do BED regulatório externo); reconstrução
ATAC com a recuperação de blocos corrompidos (código já corrigido, tracks
ainda não reconstruídas). Análise de sensibilidade já redesenhada e
re-executada sobre o código atual (`06_v1.21.2/sensitivity/sensitivity_report.md`)
- continua válida para v1.21.3 sem necessidade de repetição (opera sobre
quem passa os vetos, não sobre a distinção passes_recorded_checks/
insufficient_evidence).
