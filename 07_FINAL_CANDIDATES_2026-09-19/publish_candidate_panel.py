from pathlib import Path
import csv,json,hashlib,shutil
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'07_FINAL_CANDIDATES_2026-09-19'
assert json.loads((O/'integration_checks.json').read_text())['status']=='passed';assert json.loads((O/'sequence_validation.json').read_text())['status']=='passed'
def rows(n):return list(csv.DictReader((O/n).open(),delimiter='\t'))
panel=rows('shortlist_evidence.tsv');reserve=rows('reserve_evidence.tsv');allrows=rows('all16_decisions_and_evidence.tsv')
body='| Papel | ID | ROS_Cfam_1.0, 0-based half-open | ATAC bulk percentil | Inícios T (consenso) | RRBS cobertura / metilação observada |\n|---|---|---|---:|---:|---|\n'
for r in panel+reserve:
 role={'w01':'Prioridade 1','w11':'Alternativa condicionada','w03':'Reserva'}[r['window_id']]
 body+=f"| {role} | {r['window_id']} | {r['chrom']}:{r['start']}–{r['end']} | {float(r['atac_length_matched_percentile']):.2f} | {r['T_consensus_shared_fragment_starts']} | {float(r['rrbs_observed_fraction'])*100:.1f}% / {float(r['rrbs_mean_observed']):.2f}% |\n"
decisions='| ID | Decisão | Razão específica |\n|---|---|---|\n'+''.join(f"| {r['window_id']} | {r['decision']} | {r['reason']} |\n" for r in allrows)
note=f'''# Seleção computacional fechada para validação — Safe Harbor Canino CAR-T

**Data de fecho: 19 setembro 2026. Decisão: w01 como primeira prioridade; w11 como alternativa condicionada; w03 como reserva, sem avançar em primeira linha.** São regiões candidatas para validação. **Zero loci biologicamente validados; zero locais de corte prontos para execução.** Os intervalos de 1 kb não são sítios de inserção escolhidos ao nível da base.

Este documento consolida os resultados existentes e fecha a decisão de priorização. Não há novo score composto, relaxamento de filtros ou novo ranking automático. É uma decisão científica explícita sobre vantagens, incertezas e custo de validação. Se a definição de aprovação exigir prova de acessibilidade em CAR-T e baixa metilação, **nenhum candidato cumpre hoje esse requisito**; não é legítimo forçar um PASS.

## Painel fechado

{body}

Assembly principal: **GCF_014441545.1 ROS_Cfam_1.0**. Coordenadas acima seguem BED: início incluído, fim excluído. Para interfaces 1-based inclusive: w01 NC_051811.1:17255707–17256706; w11 NC_051836.1:27041162–27042161; w03 NC_051831.1:39713075–39714074. O crosswalk preserva assembly, orientação e limites de mapeamento. Nunca copiar números ROS para UU/CanFam6.

## 1 — w01: primeira opção para validação

Região historicamente denominada LOC119876429/LOC119872513. A janela é NC_051811.1:17255706–17256706; o gene corporal mais próximo na anotação usada é **C7H1orf21, a 110.509 bp**, também o 5′ mais próximo. O nome histórico do locus não substitui a vizinhança real. Distância corporal conferida independentemente contra o BED génico.

Vantagens: zero sobreposição com repetições anotadas no intervalo; zero hits de promotor/enhancer nas camadas EpiC interrogadas, zero overlaps de RNA não codificante/miRNA e regulação externa consultada; proxy de remapeamento **103/103 segmentos** aprovada. Correspondência UU de 1 kb, MAPQ60, NM0, sem alternativa reportada; mapeamento para CanFam6 concordante em toda a janela pelas duas rotas. Nenhuma SV sobreposta no catálogo consultado.

Limites: percentil ATAC bulk **63,75**, não acessibilidade CAR-T nem percentagem de células abertas. Dois inícios de fragmentos em dois núcleos T consenso; população T=155, um único multiome tumoral. Nas reamostragens ajustadas por profundidade, 0–2 inícios ocorrem entre percentis 2,5–97,5 dos controlos; não há fundamento para declarar enriquecimento robusto. RRBS cobre **20% da janela** e mede metilação média **60,71% nas bases observadas**; isto **não é baixa metilação demonstrada**. O máximo de 66 cães RRBS não significa 66 cães em toda a janela. Evidência do gate ATAC min/mediana/max=11/19/31 cães, sem equivalência com replicação em células T.

Variantes: **24 SNP PASS**, dos quais **10 registos com AF alternativa ≥1%**, e zero nonSNP PASS. Não tratar a janela como invariante ou escolher guia sem genótipo/seq individual. Ausência de SV no catálogo não prova ausência no animal. ARPC5 é o gene mais próximo listado no catálogo histórico de risco, a 537.787 bp; essa classificação é a do catálogo, não uma conclusão nova sobre oncogenicidade. Genes locais como C7H1orf21/TSEN15 apresentam RNA e merecem acompanhamento de perturbação.

**Decisão:** primeiro locus a testar quanto à acessibilidade no produto celular e integridade da sequência individual. Não avançar diretamente para alegação de safe harbor ou desenho final sem esses passos.

## 2 — w11: alternativa condicionada

NC_051836.1:27041161–27042161, região NPNT/TBCK. Gene corporal mais próximo **NPNT a 27.250 bp**; distância ao 5′ mais próximo 102.917 bp. Distância ao corpo e ao 5′ são perguntas diferentes e não se deve substituir uma pela outra.

ATAC bulk percentil **89,50**, dois inícios em dois núcleos T; igualmente compatível com o extremo superior das reamostragens de controlos (0–2), sem demonstração de enriquecimento T. Proxy de remapeamento 103/103 e mapeamento completo concordante. Tem **157 bp de repetição anotada**, cobertura RRBS **10%**, metilação observada **63,93%**; não é um locus demonstradamente pouco metilado. Max68 cães RRBS é máximo por local, não cobertura uniforme. Zero promotor/enhancer EpiC ou overlaps regulatórios externos nas camadas consultadas não equivale a ausência de regulação em CAR-T.

Variantes: **16 SNP PASS (8 com AF≥1%) e 2 nonSNP PASS (1 com AF≥1%)**. Sobrepõe o sinal SV de deleção catalogada com cerca de 28,9 Mb. A auditoria de genótipos anteriores encontrou heterozigosidade distribuída que questiona uma deleção constitutiva simples, mas não resolveu o sinal nem o genótipo do cão experimental. Não eliminar o aviso como artefacto sem evidência. TET2 aparece no catálogo histórico de risco a 682.054 bp; distância não prova isolamento regulatório.

**Decisão:** conservar para comparação somente com avaliação de integridade estrutural do locus no material relevante. Não o tornar primeiro apenas pelo ATAC bulk mais alto.

## Reserva — w03, região ANO2/NTF3

NC_051831.1:39713074–39714074. Bulk ATAC percentil90,95, mas **88,27% de metilação nas bases observadas**, cobertura RRBS10%, 208 bp repetitivas, um único início T e apenas998 bp de correspondência compartilhada para essa quantificação. Distância ao corpo de NTF3=25.677 bp; ao5′=97.152 bp. Não apresentou hit SV no catálogo. A prioridade inferior é justificada pelo conjunto dessas limitações; a baixa expressão de ANO2/NTF3 no RNA não supera as limitações de cromatina. Fica registada como reserva, sem recomendação de avanço imediato.

## Porque esta seleção não escolhe simplesmente o maior ATAC

Os percentis foram corrigidos para fundos de comprimento comparável. O valor ~30 que motivou a auditoria não deve continuar a ser apresentado como descrição final destes intervalos. Bulk sanguíneo, T infiltrantes num tumor e CAR-T ativadas são contextos diferentes. O objetivo não é maximizar ATAC indiscriminadamente: um sinal alto também pode coincidir com regulação endógena. Médias observadas de metilação não preenchem os gaps e zero fragmentos em dados T rasos não é veto biológico.

A evidência T usada para escolher alternativas é exploratória e foi inspecionada após geração das janelas, pelo que não representa validação independente. Das 64 alternativas Pareto anteriores, 16 receberam auditoria local adicional por apresentarem algum sinal sob consenso/união. As outras **48 permanecem guardadas, sem novo veto biológico por zero/baixo sinal T**. Esta shortlist prioriza as janelas mais documentadas, não prova que sejam as melhores possíveis no genoma nem uma classificação exaustiva de todas as64 sob os mesmos dados de variantes.

O scorer legado examinou461 regiões; a v1.21.8 manteve3 passes legados e revelou conflitos regulatórios nas regiões completas. Não se confundem esses passes com aprovação funcional. Zero overlap de uma janela não elimina os elementos regulatórios na região maior ou contactos a distância. Não se afirma aqui uma reprodução independente nova de todos os critérios do paper; a proveniência do scorer é preservada e a seleção de subintervalos é explicitamente adicional.

## Todas as16 alternativas: decisão rastreável

{decisions}

Janelas sobrepostas do mesmo locus não contam como candidatos biologicamente independentes. Por exemplo w02 sobrepõe w01 em500 bp; w10 e w12 sobrepõem w11 em750 bp. w14/w15 não são preferíveis por terem mais inícios T absolutos: têm muitos fragmentos de outros núcleos e metilação observada100% com cobertura/dadores insuficientes.

## O que falta, por ordem de decisão

| Etapa | O que permite fechar | Condição para avançar |
|---|---|---|
| Sequência e integridade no cão/material alvo | Alelos, variantes estruturais e sequência real do intervalo/flancos | Sequência resolvida; no NPNT, alerta SV resolvido antes de tratar como alvo disponível |
| Cromatina em T caninas relevantes, incluindo estado ativado/CAR-T | Acessibilidade e metilação da janela, com cobertura e replicação adequadas | Evidência reproduzível suficiente para o objetivo, sem inferir gaps como zeros; limites definidos antes de observar resultados |
| Especificidade de edição e escolha do ponto exato | Guia/ponto de corte, alelos do doador, sequência do constructo e especificidade | Desenho dependente da sequência individual e revisão de off-target; proxy de tiles não substitui isto |
| Validação pós-integração | Estrutura/número de cópias, integridade local, expressão persistente, genes vizinhos e função celular | Comparações com controlos apropriados, replicação por cão e critérios definidos antes da experiência |

Estes passos são dependentes de material/observações que não existem nos ficheiros atuais. Não podem ser honestamente marcados como concluídos por mais reprocessamento do mesmo RNA. Não foi inventado protocolo de laboratório fechado, dose, sequência de guia ou limiar arbitrário de aprovação. Se a exigência for “candidatos para levar ao laboratório”, o painel está pronto; se for “loci comprovadamente seguros”, faltam estas validações.

## Trabalho fechado nesta entrega

- Integração de16 alternativas com verificações contra tabelas de regulação, variantes, proxy de mapeabilidade e dados T originais; sem alterações de filtros.
- Dois intervalos de prioridade e uma reserva fixados, BED e crosswalk 0-based/1-based, assemblies e orientação explícitos.
- Extração independente do FASTA ROS: três sequências de1 kb sem N, hashes guardados, zero overlap de corpos génicos e distâncias recomputadas. Flancos de±1 kb fornecidos apenas como contexto, não como extensão validada da janela.
- Evidência RNA concluída: nove amostras recontadas com42.309 genes; concordância9/9 com rótulos publicados;40 genes de contexto revistos;27 testes de mates concluídos. Estes resultados não foram usados para fabricar prova de neutralidade dos loci.
- Proveniência, decisões negativas,48 alternativas não promovidas e limites de dados mantidos no pacote.

## Ficheiros para entrega

shortlist_evidence.tsv; all16_decisions_and_evidence.tsv; reserve_evidence.tsv; shortlist_ROS_0based.bed; reserve_ROS_0based.bed; coordinate_crosswalk.tsv; panel_plus_reserve_ROS_reference.fa; panel_plus_reserve_CONTEXT_flanks.fa; sequence_validation.json; integration_checks.json; input_manifest_sha256.json; remaining48_no_new_biological_veto.tsv.

As sequências são referência forward ROS, incluindo quando o mapeamento UU é reverso. Não são o genótipo do dador, não incluem guias e não selecionam uma base de corte. Estatísticas de AF são do catálogo Dog10K consultado e não previsão de compatibilidade de um cão específico. A revisão não gerou novos downloads de sequenciação.
'''
(O/'RELATORIO_FINAL_CANDIDATOS.md').write_text(note,encoding='utf-8')
brief='''# Decisão executiva — candidatos caninos CAR-T

19 setembro 2026. Seleção computacional fechada para validação experimental; zero safe harbors biologicamente validados.

1. **w01 / região LOC — primeira prioridade:** ROS_Cfam_1.0 NC_051811.1:17255706–17256706 (BED,0-based). Melhor conjunto de mapeamento, ausência de repeats anotados e cobertura RRBS entre as opções principais. ATAC bulk percentil63,75; apenas2 inícios T; metilação60,71% com20% de cobertura. Não é baixa metilação comprovada.
2. **w11 / NPNT-TBCK — alternativa condicionada:** NC_051836.1:27041161–27042161. ATAC bulk89,50;2 inícios T; metilação63,93% com10% de cobertura;157 bp repetitivos. Resolver sinal SV e genótipo antes de tratar como alvo disponível.
3. **w03 / ANO2-NTF3 — reserva:** NC_051831.1:39713074–39714074. Bulk ATAC alto não supera metilação observada88,27%, cobertura10%, repetição e sinal T escasso.

**Ainda falta para execução experimental definitiva:** sequência do cão e integridade do locus; acessibilidade/metilação no estado celular relevante; ponto de corte/guia e especificidade; validação pós-integração da expressão, estrutura, genes vizinhos e função. Nenhuma guia ou base de corte foi selecionada. Se os requisitos incluem acessibilidade CAR-T e hipometilação já demonstradas, nenhuma janela passa hoje.

O painel está pronto para decisão de validação. Evitar mais rondas de RNA como substituto dos dados em falta. Consultar RELATORIO_FINAL_CANDIDATOS.md para evidência, razões de cada alternativa e limites; BED/FASTA e hashes acompanham a entrega.
'''
(O/'LEIA_PRIMEIRO.md').write_text(brief,encoding='utf-8')
for f in ['finalize_candidate_panel.py','validate_panel_sequences.py','publish_candidate_panel.py']:shutil.copy2(f,O/f)
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - Painel final para validacao - 2026-09-19';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!important] Seleção computacional fechada —19 setembro2026\n> [[{name}]]. w01 prioridade; w11 condicionada; w03 reserva. Zero loci biologicamente validados.\n',1),encoding='utf-8')
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: seleção computacional para validação fechada\n\n[Decisão executiva](07_FINAL_CANDIDATES_2026-09-19/LEIA_PRIMEIRO.md) · [Relatório completo](07_FINAL_CANDIDATES_2026-09-19/RELATORIO_FINAL_CANDIDATOS.md).\n\nw01 prioridade1; w11 alternativa condicionada; w03 reserva. Zero safe harbors biologicamente validados; falta material/validação no contexto celular e individual para desenho definitivo.\n',encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'FINAL_COMPUTATIONAL_CANDIDATE_PANEL_20260919','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
files=[p for p in O.iterdir() if p.is_file() and p.name not in ['output_manifest_sha256.json']];(O/'output_manifest_sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
print('Decision panel, report, BED/FASTA, source audit, Obsidian and MemPalace delivered.')
