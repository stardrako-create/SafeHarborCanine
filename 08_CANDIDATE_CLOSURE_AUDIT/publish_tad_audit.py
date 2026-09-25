from pathlib import Path
import json,shutil
from mem_local import call
from audit_tad_scales import P,O
audit=json.loads((O/'tad_scale_audit.json').read_text())
lines=['\n## Atualização: sensibilidade do contexto TAD e integração de máscaras\n',
'Recalculados intervalos entre fronteiras separadamente nas escalas 100, 250 e 500 kb, além da união histórica. Em w01 e w11, nenhum dos oito intervalos contém genes do catálogo de risco usado. Este resultado é robusto às três escalas testadas, mas continua a ser uma proxy: não demonstra isolamento tridimensional em células T caninas. O código de origem descreve os dados Hi-C como destinados a scaffolding de assembly; a identidade celular relevante não fica estabelecida por esta análise.\n',
'| Janela | Escala | Intervalo ROS (0-based, end-exclusive) | Genes de risco no catálogo |\n|---|---|---|---|']
for r in audit['results']:
 lines.append(f"| {r['window']} | {r['scale']} | {r['chrom']}:{r['proxy_start']}-{r['proxy_end']} | {', '.join(r['risk_genes']) or 'nenhum'} |")
lines.extend(['\nA checklist histórica de 05_SHIP foi marcada como desatualizada. A conclusão “7/8 plenamente satisfeitos” não deve ser usada: o soft-score ATAC não valida acessibilidade em CAR-T, os intervalos entre fronteiras não validam TADs em CAR-T e já existem análises de conservação posteriores à checklist. Nenhum veto foi relaxado e as contagens anteriores foram preservadas.\n',
'Integração das duas máscaras já concluída em twenty_base_sequence_constraints.tsv: dos 981 segmentos de 20 bp por janela, 567 em w01 e 691 em w11 não têm ocorrência exata adicional nem sobreposição de REF de variantes PASS catalogadas. Isto não corresponde a guias selecionados, ausência de variantes no cão individual ou análise de mismatches/PAM.\n',
'Ainda pendente: AF específica de cada alelo normalizado; revisão adicional do SV de w11; especificidade dependente da nuclease e dados individuais. A acessibilidade e segurança funcional em CAR-T continuam por demonstrar. Nenhum worker pesado está ativo neste checkpoint; esta etapa consistiu em auditoria local concluída.\n'])
addition='\n'.join(lines)
f=O/'CURRENT.md';text=f.read_text(encoding='utf-8').replace('1. Integrar posições normalizadas das variantes com os segmentos não únicos e explicitar subintervalos limitados por esses dados, sem os promover a novos safe harbors.','1. Integração das variantes PASS com unicidade exata concluída; ver atualização abaixo. Falta refinar a AF por alelo normalizado.')
f.write_text(text+addition,encoding='utf-8')
old=P/'05_SHIP/ahmed2026_checklist_comparison.md'
banner='> CHECKLIST HISTÓRICA DESATUALIZADA — 2026-09-19: não usar a conclusão 7/8 como estado atual. ATAC soft-score e intervalos entre fronteiras Hi-C não validam acessibilidade ou isolamento em CAR-T. A conservação já foi investigada em versões posteriores. Ver ../08_CANDIDATE_CLOSURE_AUDIT/CURRENT.md.\n\n'
if not old.read_text(encoding='utf-8').startswith('> CHECKLIST HISTÓRICA'):
 shutil.copy2(old,O/'historical_ahmed_checklist_before_banner.md')
 old.write_text(banner+old.read_text(encoding='utf-8'),encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas')
dest=notes/'Safe Harbor CAR-T - Auditoria adicional de fecho - 2026-09-19.md'
dest.write_text(f.read_text(encoding='utf-8'),encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md'
if moc.exists():
 marker='> Auditoria adicional ativa em 19 setembro: painel provisório; ver [[Safe Harbor CAR-T - Auditoria adicional de fecho - 2026-09-19]]. TAD avaliado em três escalas, sem validação em CAR-T.\n\n'
 if marker not in moc.read_text(encoding='utf-8'):moc.write_text(marker+moc.read_text(encoding='utf-8'),encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'TAD_SCALE_AUDIT_20260919','content':addition,'source_file':str(dest),'added_by':'Codex'})
assert not r.get('error') and not r.get('result',{}).get('isError')
(O/'mempalace_tad_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
shutil.copy2(__file__,O/Path(__file__).name)
print('Updated CURRENT, historical checklist banner, Obsidian and MemPalace.')
