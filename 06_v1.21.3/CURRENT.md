# v1.21.3 — repeat_content/conservation preenchidos para os 3 candidatos novos

Execução em 2026-09-12. Mesma lógica de scoring de v1.21.2 (nenhuma mudança
de método); único input alterado foi o preenchimento de duas anotações que
v1.21.2 já assinalava como `missing_evidence`: `repeat_content` e
`conservation` (`max_50bp_rolling_phyloP`) para os 3 candidatos novos
(loc1 `NC_051811.1:17217706-17274511`, loc2 `NC_051820.1:22976673-23048542`,
loc3 `NC_051836.1:27013911-27071204`).

## Como foram obtidos

- **repeat_content**: RepeatMasker (Dfam 4.0, `-species dog`) sobre a FASTA
  dos 3 candidatos. loc1=29.94%, loc2=47.49%, loc3=30.23% - nenhum cruza o
  veto de 50%.
- **conservation**: liftOver ROS_Cfam_1.0→CanFam3 (`liftover/GCF_014441545.1ToCanFam3.over.chain.gz`)
  seguido de `hal2maf` sobre `241-mammalian-2020v2.hal` e `phyloP --method LRT
  --mode CONACC` contra o mesmo modelo neutro piloto de v3
  (`neutral_model_1region.mod`, ajustado a uma única região de 100kb - ver
  `05_SHIP/V3_PROGRESS_NOTES.md`; não é um modelo novo nem mais robusto).
  Máximo da média corrida de 50bp: loc1=4.01352, loc3=4.22592 - nenhum cruza
  o veto de 6.5. **loc2 não obteve conservação**: o intervalo de 71.869 bp
  não tem um bloco ortólogo único e coerente em CanFam3. Em modo `-multiple`
  o liftOver fragmenta-o em 7 pedaços com menos de 1kb cada, espalhados por
  6 cromossomas diferentes (chr15/chr9/chrX/chr14/chr6/chr5); mesmo com
  `-minMatch=0.1` (muito permissivo) o único resultado é um bloco de baixa
  confiança em chr16 que não sobrepõe nenhum dos 7 fragmentos. Isto não é
  ausência de dados - é discordância estrutural real entre as duas
  assemblies nesta região, e fica registado como `missing_evidence:
  conservation`, não forçado com um número de uma região não-sinténica.
  Note-se que loc2 também tem o repeat_content mais alto dos 3 (47.49%,
  perto do limiar de 50%) - consistente com estar numa vizinhança genómica
  mais repetitiva/instável, embora isto seja apenas circunstancial.
- Bug corrigido pelo caminho (não afeta resultados anteriores): o genoma do
  cão no HAL chama-se `Canis_lupus_familiaris`, não `dog` - o script antigo
  falhava silenciosamente neste ponto antes de qualquer phyloP correr.

## Resultado

461 candidatos: 457 excluídos, **1 com evidência insuficiente** (loc2,
falta apenas conservation agora), **3 passam as verificações registadas**
(sobe de 1 em v1.21.2):

| Candidato | Genes | Score | Rank |
|---|---|---:|---:|
| `NC_051831.1:39677324-39739751` | ANO2/NTF3 | 0.5072 | 1 |
| `NC_051836.1:27013911-27071204` | NPNT/TBCK | 0.5015 | 2 |
| `NC_051811.1:17217706-17274511` | LOC119876429/LOC119872513 | 0.3988 | 3 |

Os scores são idênticos aos já calculados em v1.21.2 (repeat_content e
conservation atuam apenas como vetos passa/falha neste scorer, não entram
na fórmula do score contínuo - confirmado empiricamente: nenhum valor
mudou). A análise de sensibilidade já corrida sobre v1.21.2
(`06_v1.21.2/sensitivity/sensitivity_report.md`) continua válida sem
necessidade de repetição: opera sobre os candidatos que passam os vetos
(hard_veto), não sobre a distinção passes_recorded_checks/
insufficient_evidence, e os 3 candidatos aqui coincidem exactamente com o
conjunto de 4 sobreviventes a p55 já reportado lá (o 4º é loc2).

**Nenhum é um safe harbor validado.** Ainda por fazer, inalterado desde
v1.21.2: reconciliar com `06_v1.21.1/CAR_T_ANALYSIS_2026-09-12.md` (PBMC vs.
células T, geometria do ponto de inserção, auditoria do BED regulatório
externo).

**Reconstrução ATAC (recuperação de blocos corrompidos) - impacto agora
verificado como negligível, prioridade reduzida**: em vez de correr a
reconstrução completa (~3h39m no cohort de 76 cães), recalculou-se
manualmente a média ponderada gate+gain da janela do atual #3
(`NC_051811.1:17217706-17274511`, que está no MESMO cromossoma do 2º caso
de corrupção encontrado em v1.20.0, `SRR27376673`/`NC_051811.1`), com os
dados reais deste cão restaurados vs. o tratamento atual de contribuição
zero. Resultado: +0,259% - na mesma ordem de grandeza do +0,238% já medido
para RIT2 em v1.20.0, novamente negligível. A reconstrução completa
continua por fazer (o código já está corrigido), mas já não é bloqueante
para confiar nos 3 candidatos atuais - fica como item de rigor/limpeza,
não como risco conhecido pendente.

Código e resultados legados (`06_v1.20.0/`, `06_v1.21.0/`, `06_v1.21.1/`,
`06_v1.21.2/`) mantidos intactos.
