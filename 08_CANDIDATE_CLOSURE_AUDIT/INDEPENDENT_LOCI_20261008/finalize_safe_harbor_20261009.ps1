$ErrorActionPreference = 'Stop'
$p = 'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino'
$a = Join-Path $p '08_CANDIDATE_CLOSURE_AUDIT'
$n = Join-Path $a 'INDEPENDENT_LOCI_20261008'
$o = Join-Path $a 'ALL_CANDIDATES_20261005'
$obs = 'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas\Safe Harbor CAR-T - ESTADO ATUAL E DECISOES.md'
$utc = [DateTime]::UtcNow.ToString('o')
$enc = New-Object System.Text.UTF8Encoding($false)
function Save-Text($path, $body) { [IO.File]::WriteAllText($path, $body, $enc) }
function Hash-File($path) { (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant() }
$report = Join-Path $o 'ANO2_REVIEW_AND_DEPENDENCIES_20261009_v7.md'
$ledger = Join-Path $n 'ANO2_closure_gate_decision_20261009_v2.json'
$manifest = Join-Path $n 'ANO2_review_delivery_manifest_20261009_v7.json'
foreach ($target in @($report,$ledger,$manifest)) { if (Test-Path -LiteralPath $target) { throw "Version already exists: $target" } }
$ncbi = Get-Content -Raw -LiteralPath (Join-Path $n 'risk_NCBI_canine_GeneID_followup_20261009_v2\receipt.json') | ConvertFrom-Json
if ($ncbi.state -ne 'completed' -or @($ncbi.still_unplaced_symbols).Count -ne 53) { throw 'Unexpected current risk receipt' }
$orth = Get-Content -Raw -LiteralPath (Join-Path $n 'ANO2_local_ortholog_risk_review_20261009\receipt.json') | ConvertFrom-Json
if ($orth.state -ne 'completed' -or @($orth.results).Count -ne 6) { throw 'Incomplete local ortholog review' }
$cov = Get-Content -Raw -LiteralPath (Join-Path $n 'risk_mapping_final_coverage_by_assembly_20261009.json') | ConvertFrom-Json
if ($cov[0].mapped -ne 2599 -or $cov[1].mapped -ne 2592) { throw 'Unexpected assembly coverage' }
$body = @'
# ANO2-39708 — resultado da etapa computacional e dependências

**Candidato adicional com suporte favorável nos testes de locus executados. Aprovação integral: NÃO emitida. Validação biológica: NÃO realizada.**

ROS_Cfam_1.0: **NC_051831.1:39708074–39709074**, coordenadas BED, zero-based, end-exclusive; 1000 bp. UU_Cfam_GSD_1.0: NC_049248.1:6741835–6742835, orientação reversa. É distinto de w01/w11, mas pertence à região da reserva w03: não contar como locus novo relativamente a w03.

## Resultados verificados

| Critério interrogado | Resultado | Cobertura e limite |
|---|---|---|
| Origem e acessibilidade | Região SHIP convergente de 62 427 bp; ATAC regional P72,70 e 20 kb P65,97, ambos ≥55 | Sangue público; não demonstra acessibilidade em células CAR-T |
| Projeção | Janela completa projetada em UU, alinhamento 1000M, NM1; correspondência CanFam6 por duas rotas | Concordância entre referências, não confirmação independente num indivíduo |
| Repetidos/especificidade da janela | 0 bp RepeatMasker dentro da janela; ROS 1000M, NM0, MAPQ60, sem alinhamentos alternativos consultados | Biblioteca, referências e modelo pesquisados |
| Conservação | 1000 bp; máximo da média de 50 bp = 1,55246 <6,5 após remoção de espécies duplicadas | 2–194 espécies não referência com ACGT/base; reanálise do mesmo alinhamento |
| RNA/TSS anotados | Sem sobreposição da janela aos genes/RNAs/exões/CDS consultados; RNA/mRNA 5′ mais próximo NTF3 a 102 152 bp ROS e 101 951 bp UU | TSS anotados, não medidos; não prova isolamento funcional |
| Regulação/cromatina | EpiC: 1000/1000 bases anotadas em 11 fontes, zero estados promotor/enhancer 1–7; sem hits TAD consultados | Proxies públicos, não ensaio de expressão ou efeito regulatório |
| SV e alelos | Sem hit no catálogo de 141 810 registos, incluindo 150 BND; cinco indels PASS normalizados; 631 segmentos de 20 bp conservados nas 32 combinações sintéticas | Não exclui SV individuais; combinações sem fase comprovada |
| Metilação | RRBS cobre 150/1000 bp; média 48,73% | Não demonstra hipometilação da janela inteira |

## Lista de risco: correspondência melhorada sem alterar a lista

Os 2654 símbolos e categorias da lista congelada foram preservados. A lista de origem usa símbolos humanos como aproximação inicial; nomes atuais/anteriores HGNC e GeneIDs caninos foram confrontados com as anotações congeladas. Nenhum símbolo foi removido para melhorar a aprovação.

| Cobertura final derivada | Símbolos colocados | Sem posição resolvida |
|---|---:|---:|
| ROS | 2599 | 55 |
| UU | 2592 | 62 |
| União ROS/UU | 2601 | 53 |

A união não corresponde à cobertura de cada genoma. Os recibos anteriores de 107 não resolvidos permanecem intactos como resultados históricos. Sete dos 53 restantes têm GeneID canino atual mas não posição nas GFF congeladas; a metadata indica outros cromossomas e foi retida sem fabricar coordenadas AR106. Os restantes 46 não tiveram correspondência atual recuperada. Uma ambiguidade HGNC foi exposta e não convertida em ortologia única.

Não foram encontrados novos genes associados à lista dentro de 300 kb. GAPDH continua o mais próximo: 879 877 bp UU; 880 859 bp ROS pela convenção de distância entre bases usada no novo recibo. O valor histórico ROS de 880 858 bp usa distância entre intervalos; diferença de convenção de 1 bp, sem alteração do gate de 300 kb. miRNA mais próximo: 1 270 499 bp UU, 1 273 355 bp ROS.

Inventário local: os mesmos seis GeneIDs nas duas referências, ANO2, NTF3 e quatro lncRNAs (LOC111092786, LOC111092811, LOC102151325, LOC102151373). A revisão NCBI por GeneID recuperou ortólogos humanos de ANO2 e NTF3, nenhum pertencente à lista congelada. Os quatro lncRNAs não têm ortólogo humano atribuído no conjunto devolvido: **correspondência de risco desconhecida**, não ausência de risco. Todas as seis respostas foram preservadas e conferidas contra contagem total/paginação. A anotação de ortologia é computacional, não validação biológica independente.

## Especificidade: duas hipóteses exploratórias

| Hipótese | Protospacer ROS BED | Sequência 5′→3′ | PAM |
|---|---|---|---|
| ANO2_TTTV23_1 | 39708151–39708174, + | TATGAGGATTTGTATCAAACACA | TTTC |
| ANO2_TTTV23_2 | 39708436–39708459, + | CCTCGGCTGTGTAATTCCTGGAG | TTTC |

Ambas tiveram zero outros hits no modelo definido em ROS completo e UU completo: 23 nt, PAMs TTTV/CTTV/TCTV/TTCV, até três substituições sem gap, ou exatamente um bulge de uma base com até duas substituições, ambas as orientações. Cada alvo pretendido foi encontrado uma vez nas coordenadas UU previstas. O índice UU corresponde aos 2197 identificadores/comprimentos da FASTA, 2 481 983 352 bases; hashes dos seis ficheiros guardados.

PAM/protospacer locais permanecem intactos nas 32 combinações sintéticas dos cinco indels PASS. Foram investigados 19 e quatro sítios de sensibilidade hipotética a uma única mudança, respetivamente: 23 projeções UU contíguas, SNPs recuperados, três indels e 40 combinações sintéticas de alelos recuperados. Nenhum criou hits no mesmo modelo. Isto não demonstra que as combinações existam em cis.

Não foram excluídos PAMs/variantes de enzima fora do modelo, múltiplos bulges, alelos fora desses 23 sítios, genomas individuais, lacunas/bases ambíguas ou cortes deslocados junto do alvo pretendido. O painel NGG anterior, com correspondências semelhantes, foi preservado. Não foi criada uma nova regra de veto por três mismatches, nem escolhida uma nuclease para obter aprovação favorável. **Estas sequências não são guias aprovados para utilização.**

## Stop real e retomada

Os workers HAL, scoring, UU e revisões desta etapa terminaram. Não existe cálculo pendente que, apenas repetido com os mesmos parâmetros, resolva a aprovação integral. A meta de pelo menos um candidato aprovado em todos os critérios continua **não atingida**.

1. Falta a **nuclease e variante realmente previstas**, para aplicar compatibilidade PAM/crRNA e o modelo de especificidade dessa enzima. Os plasmídeos históricos de Drosophila não definem o sistema canino. Retomar com essa informação, sem exigir o cão experimental para continuar a computação.
2. As posições/associações de risco não resolvidas e os quatro lncRNAs limitam a abrangência do teste. Os dados recuperados não permitem afirmar ausência de todos os riscos ou isolamento funcional. Não foram tratados como PASS nem convertidos em veto de insegurança demonstrada. Resolver esta cobertura requer anotação/associação adicional que discrimine esses casos, não novo scan idêntico.
3. Cassete, genótipo do futuro cão e ensaios de integração, expressão persistente, efeitos nos genes vizinhos e função CAR-T são dependências individuais/experimentais posteriores, separadas dos resultados computacionais.

A automação `auditar-todos-os-candidatos-safe-harbor` foi eliminada após este stop, conforme a sua regra de fecho quando não resta trabalho útil executável com os dados/parâmetros disponíveis. Não se declara o objetivo alcançado. Não foram agendadas verificações adicionais, contactados terceiros ou publicados resultados.

## Expansão e preservação

As 64 janelas alternativas pertencem a três regiões, não a 64 loci independentes; as 461 regiões originais e exclusões históricas foram preservadas. A conservação dos independentes 1059/1266/1900 ficou concluída: 0,50206 / 2,21532 / 0,46692, 1000 bp cada. Os seus gaps estão fora da geometria SHIP original de 50–75 kb; 1059 também falha ATAC do gap e 1059/1900 não são convergentes. Não foram usadas extensões favoráveis escolhidas depois de ver os resultados para substituir esse gate. w11 continua fora do avanço por estrutura não resolvida, sem alegar insegurança demonstrada.

CURRENT, decisões, ficha ANO2 e Obsidian atualizados. MemPalace voltou a falhar com MCP SSE 404 no fecho; fila local mantida, **sem alegar sincronização**. Scripts, tentativas falhadas, alinhamentos brutos e versões anteriores preservados.

## Recibos e fontes

Em `../INDEPENDENT_LOCI_20261008/`, o manifesto `ANO2_review_delivery_manifest_20261009_v7.json` guarda SHA256 dos recibos e deste relatório. Fontes/resultados principais:

- `ANO2_two_guides_full_UU_20261009/results.json` e `full_UU_index_integrity_20261009_v2/receipt.json`;
- `ANO2_TTTV23_offsite_UU_and_catalogue_20261009/results.json`, `ANO2_offsite_three_indels_20261009/receipt.json`, `ANO2_offsite_catalogue_combinations_20261009/receipt.json`;
- `ANO2_UU_annotation_crosscheck_20261009/receipt.json` (GFF oficial, MD5 e SHA256);
- `risk_HGNC_rename_audit_20261009/receipt.json`, `risk_NCBI_canine_GeneID_followup_20261009_v2/receipt.json`, `risk_mapping_final_coverage_by_assembly_20261009.json`;
- `ANO2_local_gene_inventory_20261009.json`, `ANO2_local_ortholog_risk_review_20261009/receipt.json`;
- `remote_conservation_taxa_review_20261009/receipt.json` e recibos de conservação ANO2 referenciados nas versões anteriores.

Anotação: [NCBI AR106](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Canis_lupus_familiaris/106/). Nomenclatura: [campos HGNC](https://www.genenames.org/help/statistics-and-downloads/). GeneIDs: [NCBI GeneInfo Mammalia](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/). Ortólogos: [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/tutorials/download-ortholog-dataset/). Snapshots, respostas e hashes locais são a proveniência das contagens acima.
'@
Save-Text $report ($body + "`n`nFecho registado em UTC: $utc`n")
$sources = @(
 'ANO2_two_guides_full_UU_20261009\results.json',
 'full_UU_index_integrity_20261009_v2\receipt.json',
 'ANO2_UU_annotation_crosscheck_20261009\receipt.json',
 'ANO2_offsite_catalogue_combinations_20261009\receipt.json',
 'risk_HGNC_rename_audit_20261009\receipt.json',
 'risk_NCBI_canine_GeneID_followup_20261009_v2\receipt.json',
 'risk_mapping_final_coverage_by_assembly_20261009.json',
 'ANO2_local_gene_inventory_20261009.json',
 'ANO2_local_ortholog_risk_review_20261009\receipt.json',
 'remote_conservation_taxa_review_20261009\receipt.json'
) | ForEach-Object { $q = Join-Path $n $_; [ordered]@{path=$q;sha256=(Hash-File $q)} }
$decision = [ordered]@{
 utc=$utc; candidate='ANO2-39708'; locus_filters_supported_in_interrogated_data=$true;
 complete_audit_approval=$false; biological_validation=$false; newly_fully_approved_count=0;
 actual_nuclease_selected=$false; blocking_input='Actual nuclease and variant for enzyme-specific guide acceptance';
 risk_mapping_union=@{mapped=2601;unresolved=53;total=2654}; risk_mapping_by_assembly=$cov;
 local_lncRNA_risk_correspondence_unknown_count=4; unknowns_are_not_pass_or_demonstrated_veto=$true;
 workers_this_batch_completed=$true; pending_worker_count=0;
 stop_reason='Defined reference scans and bounded annotation follow-ups complete; actual nuclease missing and remaining annotation/function coverage not resolved by these data';
 automation=@{id='auditar-todos-os-candidatos-safe-harbor';deleteStatus='deleted';goal_achieved=$false};
 mempalace=@{sync=$false;end_retry='MCP SSE HTTP404';local_queue_preserved=$true};
 source_hashes=$sources;report=@{path=$report;sha256=(Hash-File $report)}
}
Save-Text $ledger ($decision | ConvertTo-Json -Depth 14)
$front = @"
## Estado operacional — 09/10/2026, etapa concluída com stop documentado
Este resumo substitui os estados operacionais históricos abaixo. **Zero novos candidatos integralmente aprovados; nenhuma validação biológica.**

- ANO2-39708: ROS NC_051831.1:39708074–39709074 (BED), candidato adicional com suporte favorável nos testes executados. Distinto de w01/w11, mesma região da reserva w03. Relatório atual: ALL_CANDIDATES_20261005/ANO2_REVIEW_AND_DEPENDENCIES_20261009_v7.md; ledger v2 e manifesto v7 em INDEPENDENT_LOCI_20261008.
- HAL, scoring, scans completos ROS/UU e revisões desta etapa concluídos. Duas hipóteses TTTV23 com zero outros hits no modelo explícito; não são guias aprovados. Índice UU integral conferido. Não duplicar trabalhos terminados.
- Lista congelada 2654 preservada: união ROS/UU 2601 colocados, 53 não resolvidos; por genoma ROS 2599/55 e UU 2592/62. Não houve novo hit de risco dentro de 300 kb. ANO2/NTF3 com ortólogos humanos fora da lista; quatro lncRNAs com correspondência de risco desconhecida. Ausência de correspondência não é PASS nem insegurança demonstrada.
- Stop real: nuclease/variante não indicada para fechar especificidade da intervenção; lacunas de anotação/associação não permitem afirmar ausência de todos os riscos. Retomar com parâmetro/dado novo que resolva estas questões; não exigir o cão experimental para continuar a computação.
- Conservação 1059/1266/1900 concluída (0,50206/2,21532/0,46692); geometria SHIP original não satisfeita. Não relaxar critérios para atingir a meta. w11 fora do avanço por estrutura não resolvida; w01 sem validação funcional.
- Automação de 30 minutos auditar-todos-os-candidatos-safe-harbor eliminada no stop; objetivo não atingido. Nenhuma verificação adicional agendada. Dependências individuais/experimentais descritas no relatório.
- CURRENT e Obsidian atualizados; MemPalace novamente indisponível (MCP SSE 404), fila local preservada. Sem contactos externos ou publicação. Fecho UTC: $utc.

"@
$note = @"

## Fecho documentado — $utc
Relatório ANO2 v7 e decisão v2 concluídos. Cobertura de risco derivada: união 2601/2654, 53 não resolvidos; ROS 2599/55, UU 2592/62. Seis genes locais revistos por GeneID/ortólogos: ANO2/NTF3 sem associação à lista congelada, quatro lncRNAs sem correspondência humana atribuída, mantendo risco desconhecido. Zero novas aprovações integrais. Workers desta etapa terminados. Stop: nuclease/variante real em falta e cobertura de anotação/associação insuficiente para afirmações universais; não é dependência do cão para executar os scans. Automação auditar-todos-os-candidatos-safe-harbor eliminada conforme regra de fecho, sem alegar meta alcançada. MemPalace reavaliado no fecho: SSE404; fila local preservada. Inputs, outputs e contagens históricas intactos; nenhum contacto externo/publicação.
"@
$stateFiles = @((Join-Path $p 'CURRENT.md'),(Join-Path $a 'CURRENT_STATE_AND_DECISIONS.md'),$obs)
foreach ($f in $stateFiles) {
 $old = [IO.File]::ReadAllText($f)
 $marker = '### Histórico preservado a partir daqui'
 $idx = $old.IndexOf($marker)
 $begin = $old.IndexOf('## Estado operacional')
 if ($begin -lt 0 -or $idx -lt $begin) { throw "Missing bounded front markers in $f" }
 $prefix = $old.Substring(0,$begin)
 Save-Text $f ($prefix + $front + $old.Substring($idx) + $note + "`n")
}
foreach ($f in @((Join-Path $a 'CANDIDATE_EXPANSION_20261005.md'),(Join-Path $o 'CANDIDATO_ANO2_39708.md'),(Join-Path $a 'MEMPALACE_SYNC_PENDING_20260922.md'))) {
 [IO.File]::AppendAllText($f,$note + "`n",$enc)
}
$manifestBody = [ordered]@{utc=$utc;complete_audit_approval=$false;biological_validation=$false;source_hashes=$sources;report=@{path=$report;sha256=(Hash-File $report)};decision=@{path=$ledger;sha256=(Hash-File $ledger)};state_files=@($stateFiles | ForEach-Object { @{path=$_;sha256=(Hash-File $_)} });automation_deleted=$true;mempalace_synced=$false}
Save-Text $manifest ($manifestBody | ConvertTo-Json -Depth 14)
foreach ($f in @($report,$ledger,$manifest)) { if ((Get-Item -LiteralPath $f).Length -eq 0) { throw "Empty output $f" } }
@{report=$report;ledger=$ledger;manifest=$manifest;utc=$utc;approved=0;unresolved_union=53} | ConvertTo-Json
