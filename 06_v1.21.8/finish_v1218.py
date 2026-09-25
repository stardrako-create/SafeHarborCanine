from pathlib import Path
import json,hashlib,shutil,datetime
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.21.8';stamp=datetime.datetime.now().isoformat(timespec='seconds')
summary=json.loads((O/'scoring_comparison_summary.json').read_text());local=json.loads((O/'local_refresh_summary.json').read_text());assert json.loads((O/'rebuild_targeted_validation.json').read_text())['status']=='passed';assert json.loads((O/'rebuild_global_validation.json').read_text())['status']=='passed';assert json.loads((O/'scoring_execution.json').read_text())['returncode']==0
test=P/'scripts/tests/test_regression_v1218.py';s=test.read_text();s=s.replace("if __name__=='__main__':unittest.main()\n",'');s+="\nif __name__=='__main__':unittest.main()\n";test.write_text(s)
(O/'tests').mkdir(exist_ok=True)
for p in (P/'scripts/tests').glob('test_*.py'):shutil.copy2(p,O/'tests'/p.name)
for name in ['build_mother_track_v2.py','score_ship_candidates_v2.py','bw_utils.py']:shutil.copy2(P/'scripts'/name,O/'code'/name)
for p in Path('.').glob('*v1218.py'):shutil.copy2(p,O/p.name)
(O/'test_summary.json').write_text(json.dumps(dict(status='passed',tests=18,command="python -m unittest discover -s scripts/tests -v",observed_at=stamp,coverage=['RRBS missingness','matching window statistics','gene 5-prime distances','missing weighted score components','unrecoverable read raises','recoverable read retains data','terminal partial-bin denominator','written mean/variability gaps for zero/one/two dogs']),indent=2))
note='''# v1.21.8 — reconstrução ATAC e reavaliação verificadas

Fecho desta etapa técnica em STAMP. Não é fecho da validação biológica do projeto.

## Resultado

A reconstrução global e o rescoring foram concluídos e revistos. **As mesmas três regiões passam o scorer legado; nenhuma é um safe harbor CAR-T validado.** Na avaliação integrada, as três regiões completas têm anotações regulatórias que exigem resolução local. A janela específica de LOC não tem sobreposição EpiC, mas a região que a contém tem.

| Região | Score legado atualizado | Percentil ATAC da região | Percentil da janela local anterior de 1 kb | Promotor/enhancer EpiC na região completa |
|---|---:|---:|---:|---:|
| ANO2/NTF3 | 0,5073 | 72,70 | 93,98 | 9.616 bp |
| LOC119876429/LOC119872513 | 0,3984 | 57,43 | 96,03 | 5.349 bp |
| NPNT/TBCK | 0,5015 | 57,63 | 96,57 | 14.800 bp |

Percentis regionais e locais respondem a perguntas diferentes, com os respetivos fundos de comprimento correspondente. Janelas escolhidas por máximos anteriores continuam exploratórias. Estes números não medem percentagem de cromatina aberta e não validam células T/CAR-T.

## Verificação da reconstrução

76 cães, 376 sequências e 2.396.858.295 bp, em bins de 25 bp. Foram verificadas as invariantes de **95.874.515 bins** nos caches e as coberturas dos BigWigs finais. Fundo de ganho reproduzido: 0,08076513558626175. Média observada em 2.328.830.335 bp e variabilidade com ≥2 cães em 2.312.901.198 bp. Ausência de cães com evidência mantém-se gap, não zero medido.

Quatro blocos foram recuperados de raw BigWigs. A conversão raw→CPM foi conferida em múltiplas sub-regiões legíveis de cada bloco; as sub-regiões ainda ilegíveis no CPM foram registadas e não usadas como falsa confirmação. Em 20 intervalos independentes (815 bins), o cálculo direto por base de média, IQR e número de cães concorda com os novos tracks, incluindo as três janelas e regiões de recuperação. Isto verifica os cálculos e escalas, não certifica os alinhamentos/experimentos originais.

**18 testes de regressão passaram**, incluindo um teste que escreve e volta a ler BigWigs com zero, um e dois cães: média ausente com zero, variabilidade ausente com menos de dois, e valores corretos quando há evidência. Corrigidos no código principal: abortar perante bloco irrecuperável em vez de usar zero; explicitar componentes de score ausentes; excluir bases de padding da média do último bin; preservar gaps na escrita de médias. O construtor rápido usado nesta execução e o código principal estão preservados. Os seus caminhos de leitura diferem; a equivalência numérica foi verificada nos intervalos de teste, não por reexecutar o construtor lento integralmente.

## Comparação com a versão anterior

461 regiões: 457 excluídas, uma com evidência insuficiente e três passes legados. **Zero mudanças no estado final ou no hard veto agregado.** Houve mudança em 205 scores finais, 246 componentes de estabilidade ATAC e 13 decisões do filtro isolado de acessibilidade; essas 13 não alteraram o estado agregado porque outros critérios permanecem ativos. As diferenças completas estão em scoring_changes.tsv. A nova regra de componentes ausentes não alterou o estado final destes 461 registos, mas fecha um comportamento incorreto para entradas incompletas.

Foi criado candidates_integrated_review_v1218.tsv: conserva o resultado legado e acrescenta a revisão regulatória regional e ausência de validação funcional. Não transforma um PASS de score em aprovação de safe harbor. Os conflitos regionais não eliminam automaticamente todos os subintervalos possíveis.

## Alternativas locais

As 2.091 janelas foram atualizadas com os novos tracks, cobertura e contagem de cães. As três janelas selecionadas anteriormente conservam cobertura ATAC observada de 100%. Cães por base/bin, mínimo/mediana/máximo: ANO2 26/33/42; LOC 16/34/57; NPNT 12/28/62. Isto é presença acima do gate, não independência entre bins nem medida validada de atividade em T.

Para orientar trabalho futuro, há 306 janelas de 1 kb sem sobreposição regulatória registada no filtro EpiC + BED externo e com covariáveis disponíveis para comparação: 116 ANO2, 115 LOC, 75 NPNT. Um conjunto de 64 janelas apresenta compromissos não dominados dentro de cada região entre maior ATAC, maior cobertura RRBS, menos repetições e menor metilação observada. Nenhum peso/limiar novo foi escolhido para forçar uma passagem. Não são 64 loci independentes: muitas janelas sobrepõem-se. A maior parte não tem as verificações de variantes, remapeamento e dados T feitas nas três janelas anteriores. Não são nomeações finais para inserção. Ver local_tradeoff_frontier_exploratory.tsv.

## O que permanece aberto

1. Identidade celular/QC completos no multiome, replicação em T caninas saudáveis/ativadas e contexto CAR-T. Os resultados anteriores de 0/1/0 fragmentos no grupo restrito de 58 núcleos permanecem insuficientes; a reconstrução do sangue não altera esses dados independentes.
2. Calibração adequada da conservação e explicação completa das pequenas diferenças de reprodução. O teste de escala em segunda região é sensibilidade, não validação de um modelo neutro.
3. Confirmação da grande SV NPNT e sequência individual; os catálogos populacionais não substituem genotipagem do animal relevante.
4. Resolução dos conflitos regulatórios e avaliação completa das alternativas antes de escolher qualquer ponto de inserção.
5. Demonstração experimental de expressão estável e neutralidade funcional após integração. Não está resolvida por este pipeline.

Por isso, a etapa de reconstrução/validação/rescoring está concluída; o objetivo de demonstrar um safe harbor canino para CAR-T não está. A conclusão continua a ser shortlist exploratória com riscos explícitos, sem vencedor validado.

## Ficheiros auditáveis

ATAC/build_manifest.json; rebuild_targeted_validation.json; rebuild_global_validation.json; test_summary.json; scoring_command.json; scoring_execution.json; scoring.log; scoring_changes.tsv; scoring_comparison_summary.json; candidates_scored_v1218.tsv; candidates_integrated_review_v1218.tsv; local_windows_evidence_status.tsv; local_refresh_summary.json; local_tradeoff_frontier_exploratory.tsv. O relatório da reunião original fica preservado como retrato anterior; o relatório atualizado contém o estado agora verificado.
'''.replace('STAMP',stamp)
(O/'CURRENT.md').write_text(note,encoding='utf-8')
report=(O/'RELATORIO_REUNIAO_2026-09-14.md').read_text(encoding='utf-8')
# Preserve the meeting snapshot; make an updated full report with an explicit current status.
report=report.replace('**Estado:** resultados auditados até v1.21.7; reconstrução ATAC v1.21.8 concluída, ainda sem validação final e novo scoring.','**Estado atualizado:** reconstrução, validação numérica e novo scoring v1.21.8 concluídos; validação biológica ainda pendente.')
report=report.replace('A última comparação completa anterior à reconstrução atual avaliou','A comparação completa atual avaliou').replace('Estes números referem-se ao scorer legado; não são a contagem de safe harbors validados nem o resultado de um novo scoring v1.21.8.','O rescoring v1.21.8 preservou estes estados. São resultados dos critérios legados, não contagem de safe harbors validados.')
report=report.replace('| 94,1 |','| 93,98 |').replace('| 96,07 |','| 96,03 |').replace('| 96,63 |','| 96,57 |')
report=report.replace('Os percentis são dos cálculos locais anteriores, contra 3.000 janelas de fundo com comprimento correspondente, e ainda não foram recalculados com os novos tracks.','Os percentis foram recalculados com os novos tracks, usando 3.000 tentativas de janelas de fundo por comprimento e excluindo janelas sem média observada. Os intervalos candidatos são os selecionados anteriormente.')
a=report.index('## 6. Atualização de hoje: v1.21.8');b=report.index('## 7. O que falta',a)
section=note[note.index('## Verificação da reconstrução'):note.index('## O que permanece aberto')]
report=report[:a]+'## 6. Atualização verificada: v1.21.8\n\n'+section+'\n'+report[b:]
report=report.replace('| Imediata, computacional | Validar v1.21.8, recuperações e alterações; recalcular scoring e rever diferenças | Uma versão de resultados consistente, rastreável e revista |','| Concluída nesta etapa | Validação numérica da v1.21.8, testes específicos, rescoring e comparação | 18 testes e verificações de recuperação/caches passaram; shortlist legada preservada |')
report=report.replace('- «A reconstrução terminada valida os resultados finais» — a execução terminou, mas validação e novo scoring estão pendentes.','- «A validação numérica prova segurança biológica» — a reconstrução e o rescoring foram verificados, mas acessibilidade/neutralidade em CAR-T continuam por demonstrar.')
updated=O/'RELATORIO_ATUALIZADO_2026-09-14.md';updated.write_text(report,encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.21.8 validada e reavaliada - 2026-09-14';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');(notes/'Safe Harbor CAR-T - Relatorio atualizado - 2026-09-14.md').write_text(report,encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC_before_validation.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!important] v1.21.8 — reconstrução e rescoring verificados\n> [[{name}]]; [[Safe Harbor CAR-T - Relatorio atualizado - 2026-09-14]]. As mesmas três regiões legadas; conflitos regulatórios e validação CAR-T pendentes.\n',1),encoding='utf-8')
shutil.copy2(P/'CURRENT.md',O/'previous_CURRENT_before_validation.md');(P/'CURRENT.md').write_text('# Estado atual: v1.21.8 verificada\n\n[Checkpoint](06_v1.21.8/CURRENT.md) · [Relatório atualizado](06_v1.21.8/RELATORIO_ATUALIZADO_2026-09-14.md).\n\nReconstrução ATAC, testes, rescoring e comparação concluídos. Mesmas três regiões legadas; nenhuma validação de safe harbor CAR-T. Questões biológicas e de anotação permanecem explícitas.\n',encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.21.8_Validated_Rescore','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'validation_mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
hashes={}
for p in list(O.glob('*.py'))+list((O/'code').glob('*.py'))+list((O/'ATAC').glob('*.bw'))+list(O.glob('*.tsv')):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
 hashes[str(p.relative_to(O))]=h.hexdigest()
(O/'validated_output_hashes.json').write_text(json.dumps(hashes,indent=2));print('Validated stage and updated report saved; overall biological validation still open.')
