from pathlib import Path
import csv,json,shutil,hashlib
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.0'
rows=list(csv.DictReader((O/'frontier_T_evidence_v1220.tsv').open(),delimiter='\t'));summary=json.loads((O/'frontier_T_summary.json').read_text());assert len(rows)==64 and all(x['HTSlib_match'] for x in summary['validation'])
for r in rows:
 assert 0<float(r['source_mapping_agreement_fraction'])<=1
 assert int(r['T_consensus_shared_fragment_starts'])<=int(r['T_union_shared_fragment_starts'])
 if r['canfam6_mapping_status']=='two_routes_unique_identical_contiguous':
  assert int(r['shared_mapped_bp'])==1000
  for mask in ['T_consensus','T_union','all_QC']:assert int(r[mask+'_fragment_starts'])==int(r[mask+'_shared_fragment_starts'])
positive=[r for r in rows if int(r['T_consensus_shared_fragment_starts'])>0];assert len(positive)==13
with (O/'positive_windows_exploratory.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(positive)
note='''# Safe Harbor CAR-T — v1.22.0: ATAC nas alternativas locais

14 de setembro de 2026. **A análise foi alargada às 64 janelas alternativas. Há sinais T esparsos fora das três janelas anteriores, mas continuam a faltar suporte replicado e validação funcional. Não há vencedor demonstrado.** O scoring global continua a ser o da v1.21.8; esta versão acrescenta evidência local.

## O que foi feito

Mapeámos as 64 janelas de 1 kb por duas vias: ROS→CanFam3→CanFam6 e inversão da cadeia CanFam6→ROS. Exigimos concordância por base e unicidade do destino. Em 37 janelas, o mapeamento é integral e contínuo. Nas outras 27, há pequenas lacunas, diferenças entre vias ou descontinuidades: contamos separadamente apenas inícios em bases de destino com concordância, guardando explicitamente a cobertura. Não imputámos zero às bases não avaliáveis.

Foram consultados três intervalos regionais de fragmentos do [multiome público GSE244116](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244116), reutilizando o cache. As consultas foram bracketed por registos ordenados antes/depois dos intervalos. A descompressão independente com HTSlib reproduziu os 2.027, 4.613 e 5.167 registos guardados para os intervalos de LOC, ANO2 e NPNT. A completude das contagens de início depende da ordenação do ficheiro; não é uma contagem exaustiva de fragmentos longos sobrepostos cujo início esteja fora da consulta.

Usámos os mesmos grupos definidos por RNA na v1.21.9: 155 núcleos na interseção das quatro análises, 157 na união, 5.849 núcleos QC no total. Os resultados por janela são idênticos entre interseção e união. São núcleos de um único tumor, não réplicas biológicas ou T caninas saudáveis validadas.

## Resultado

| Região | Alternativas | Mapeamento integral/contínuo | Menor cobertura de concordância nas alternativas | Janelas com ≥1 início no grupo T, contando as bases concordantes | Fragmentos distintos / núcleos distintos em todas as alternativas da região |
|---|---:|---:|---:|---:|---:|
| ANO2/NTF3 | 22 | 13 | 99,7% | 1 | 1 / 1 |
| LOC119876429/LOC119872513 | 12 | 5 | 92,3% | 3 | 3 / 3 |
| NPNT/TBCK | 30 | 19 | 93,6% | 9 | 9 / 9 |

As 13 janelas com algum sinal não são 13 descobertas independentes: há sobreposição e reutilização dos mesmos fragmentos entre janelas. Os totais da última coluna removem essa duplicação dentro de cada região. As bases sem correspondência continuam sem avaliação; contagens nas bases concordantes não equivalem necessariamente a uma janela integral de 1 kb.

Entre as 37 janelas integralmente mapeadas, têm algum sinal zero de 13 em ANO2, duas de cinco em LOC e oito de 19 em NPNT. A inclusão cuidadosa das bases concordantes das restantes 27 acrescenta uma janela positiva em cada região.

## Janelas úteis para decidir a próxima verificação

Estas são exemplos de compromissos, não locais aprovados para edição. Coordenadas ROS, 0-based half-open. Todas pertencem à frente exploratória anteriormente definida sem usar estas contagens T; priorizá-las agora com estes mesmos dados continua a exigir validação independente.

| Janela | Inícios no grupo T de 155 | Inícios em todos os 5.849 QC | Percentil ATAC no sangue | Fração RRBS observada | Máximo da média phyloP publicada em 50 bp |
|---|---:|---:|---:|---:|---:|
| LOC: NC_051811.1:17255706–17256706 | 2 | 18 | 63,75 | 20,0% | 0,342 |
| NPNT: NC_051836.1:27041161–27042161 | 2 | 24 | 89,50 | 10,0% | 0,314 |
| NPNT: NC_051836.1:27050411–27051411 | 4 | 490 | 66,61 | 3,9% | 1,165 |

Estas três janelas têm mapeamento integral concordante. A janela NPNT com quatro inícios T tem também 490 inícios no conjunto QC: quatro não demonstra enriquecimento específico em T. Não foi feita normalização pela profundidade genómica de fragmentos nem inferência estatística de enriquecimento. A cobertura RRBS de 3,9% é particularmente limitada; o resto não pode ser chamado desmetilado.

A primeira janela LOC oferece um ponto concreto para aprofundar a revisão local pela menor conservação e ausência de sobreposição regulatória no filtro que definiu a frente. Duas moléculas em dois núcleos continuam a ser suporte muito escasso. A janela NPNT com percentil sanguíneo 89,50 permite estudar outro compromisso, também com só dois núcleos e RRBS incompleto. Nenhuma está sujeita ainda a todas as verificações individuais de variantes/mapeabilidade efetuadas nas três janelas originais. A chamada SV regional NPNT mantém a interpretação revista na v1.21.9: catalogada e questionada por heterozigotia distribuída, ainda por resolver.

## O que isto muda

Os zeros das três janelas anteriores não se generalizam a todas as alternativas. Há sinais focais noutros subintervalos, mas a dimensão da evidência continua insuficiente para afirmar acessibilidade robusta em T ou segurança para CAR-T. Procurar apenas o máximo de ATAC no sangue não seleciona necessariamente a janela com melhor suporte na população T exploratória.

A próxima etapa computacional útil é aplicar às alternativas priorizadas a revisão de variantes locais, mapeabilidade e contexto regulatório com a mesma exigência das janelas originais, e completar o QC/identidade do multiome. Para fechar a conclusão biológica, continuam indispensáveis dados de T caninas relevantes com replicação e validação da integração/neutralidade funcional. Não se deve baixar filtros para promover estes sinais escassos a candidatos validados.

## Verificação e ficheiros

Foram verificadas 64 linhas, limites de cobertura, identidade das contagens interseção/união e equivalência das contagens por bases concordantes com as contagens integrais nas 37 janelas aplicáveis. A descompressão HTSlib coincide com os registos consultados. Os critérios anteriores de seleção e scores foram preservados.

`frontier_T_evidence_v1220.tsv` contém todas as 64 janelas, cobertura e contagens integrais/parciais separadas; `positive_windows_exploratory.tsv`, as 13 com algum sinal; `frontier_shared_mapping_blocks.json`, os blocos auditáveis; `frontier_mapping_summary.json`, o critério de mapeamento integral; `frontier_T_summary.json`, a análise que acrescenta as bases concordantes das restantes janelas; `regional_fragment_queries.json`, a proveniência das consultas. Os intervalos com concordância parcial não receberam aprovação de mapeamento integral.
'''
differences=[r for r in rows if r['T_consensus_shared_fragment_starts']!=r['T_union_shared_fragment_starts']];assert len(differences)==3 and all(r['genes']=='ANO2/NTF3' for r in differences)
note=note.replace('Os resultados por janela são idênticos entre interseção e união.','Em 61 das 64 janelas, as contagens são iguais entre interseção e união. Em três janelas sobrepostas ANO2 (inícios ROS 39729824, 39730074 e 39730324), a união acrescenta um início por janela onde a interseção tem zero. Assim, são 13 janelas positivas na interseção e 16 na união; a tabela abaixo refere-se à interseção de 155.')
note=note.replace('identidade das contagens interseção/união','inclusão das contagens da interseção nas da união e identificação explícita das três diferenças')
(O/'CURRENT.md').write_text(note,encoding='utf-8')
report=note+'\n\n---\n\n# Contexto anterior preservado — relatório v1.21.9\n\nA análise acima atualiza a acessibilidade das alternativas, que neste checkpoint anterior ainda estava pendente.\n\n'+(P/'06_v1.21.9/RELATORIO_ATUALIZADO_2026-09-14.md').read_text(encoding='utf-8')
(O/'RELATORIO_ATUALIZADO_2026-09-14.md').write_text(report,encoding='utf-8')
for name in ['map_frontier_v1220.py','prepare_queries_v1220.py','quantify_frontier_v1220.py','finish_v1220.py']:shutil.copy2(Path(name),O/name)
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.0 ATAC nas alternativas - 2026-09-14';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8')
shared=notes/'Safe Harbor CAR-T - Relatorio atualizado - 2026-09-14.md';shutil.copy2(shared,O/'previous_Obsidian_report.md');shared.write_text(report,encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!important] Estado mais recente: v1.22.0\n> [[{name}]] — 64 alternativas interrogadas nas bases concordantes; 37 integralmente mapeadas; sinais T escassos em 13 janelas sobrepostas, sem vencedor validado. Entradas seguintes são histórico.\n',1),encoding='utf-8')
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.22.0\n\n[Relatório atualizado](06_v1.22.0/RELATORIO_ATUALIZADO_2026-09-14.md) · [Evidência das alternativas](06_v1.22.0/frontier_T_evidence_v1220.tsv).\n\n64 alternativas avaliadas nas bases concordantes; 37 integralmente mapeadas. Sinais T escassos, sem candidato validado. Scoring global permanece v1.21.8.\n',encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.0_Alternative_Windows_T_Evidence','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'validation_summary.json').write_text(json.dumps(dict(status='passed',windows=64,full_contiguous=37,positive_shared_windows=13,positive_union_shared_windows=16,consensus_union_differing_windows=3,HTSlib_queries_verified=3,checks=['coverage bounds','consensus counts bounded by union; three differing windows identified','full/shared count equality on all37 full mappings','64 distinct original windows preserved']),indent=2))
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='output_hashes.json'};(O/'output_hashes.json').write_text(json.dumps(hashes,indent=2))
assert dest.read_text(encoding='utf-8')==note and shared.read_text(encoding='utf-8')==report
print('v1.22.0 saved and checked: project, Obsidian, MemPalace, report. No validated safe harbor.')
