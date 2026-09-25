from pathlib import Path
import json,shutil
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.23.0';s=json.loads((O/'protocol_structure_audit.json').read_text());assert len(s['metadata'])==9 and len(s['read_structure'])==18
r2=[x for x in s['read_structure'] if x['mate']==2];n=sum(x['n'] for x in s['read_structure']);assert n==1800000
missing_submitted=all(not r.get('submitted_ftp') and not r.get('submitted_format') for x in s['metadata'] for r in x['submitted'])
lengths=sorted({k for x in s['read_structure'] for k in x['length_histogram']});phreds=sorted({k for x in s['read_structure'] for k in x['phred_histogram']})
note=f'''# v1.23.0 — auditoria de protocolo e estrutura dos reads

19 setembro 2026. **A descrição do kit e a estrutura dos FASTQs disponíveis não estão reconciliadas.** Esta etapa delimita a dúvida com dados; não identifica qual kit foi efetivamente usado nem invalida automaticamente a contagem anterior.

## Fontes verificadas

Descarregados os XML dos nove experimentos ENA ligados aos runs do manifesto. Todos declaram RNA-Seq, layout PAIRED e SMART-Seq v4 3’ DE. Essa concordância documental pode refletir metadados copiados do mesmo depósito; não é confirmação independente do procedimento experimental. Exemplo: [SRX22428478](https://www.ebi.ac.uk/ena/browser/api/xml/SRX22428478). A descrição concorda com o [GEO GSM7887650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7887650).

O manual original da Takara, versão 092618, páginas 8 e 21, descreve R2 com um índice inline de seis bases seguido de poly(T), destinado a desmultiplexagem. Fonte primária de autoria do fabricante, consultada numa [cópia do manual](https://device.report/m/a1cff4f8333dfe571d210a1f84e734b6cdb0f611dc09ffa2d1e02c4b1109f9f6.pdf); a ligação atual do fabricante não ficou acessível nesta sessão. Esta descrição não deve ser extrapolada a todos os protocolos SMART-Seq nem aplicada cegamente a dados já processados.

## Verificação dos ficheiros disponíveis

Inspecionados os primeiros 100.000 reads de cada um dos 18 FASTQs ENA: {n:,} reads, amostra inicial determinística, não amostragem aleatória do ficheiro inteiro. Comprimentos encontrados: {', '.join(lengths)} bases. Valores Phred observados: {', '.join(phreds)}.

Em R2, a fração com prefixo exatamente igual a um dos 12 índices do manual foi {min(x['known_inline_prefix_fraction'] for x in r2)*100:.3f}–{max(x['known_inline_prefix_fraction'] for x in r2)*100:.3f}%. A fração com pelo menos 16 T nas posições 7–26 foi {min(x['polyT_80pct_positions7_26_fraction'] for x in r2)*100:.3f}–{max(x['polyT_80pct_positions7_26_fraction'] for x in r2)*100:.3f}%. A assinatura conjunta índice+poly(T) ocorreu em {min(x['joint_index_polyT_fraction'] for x in r2)*100:.3f}–{max(x['joint_index_polyT_fraction'] for x in r2)*100:.3f}% dos reads R2. O limiar 16/20 é uma definição exploratória explícita, não uma especificação do fabricante nem teste formal do kit.

Não se observa a assinatura esperada como estrutura dominante no início destes reads. Isto é compatível com processamento prévio, uma preparação diferente ou uma descrição incompleta; os dados não distinguem essas hipóteses. Não prova qual delas ocorreu. As contagens GeneCounts e o mapeamento elevado são evidência operacional adicional, mas não identificam o protocolo nem confirmam por si o papel de cada mate.

Consultados também campos submitted_format/ftp/md5/bytes da ENA para os nove runs. Ausentes os campos de formato e ligação de ficheiros originalmente submetidos em todas as respostas: {missing_submitted}. Campos vazios não demonstram ausência dos originais nem estabelecem proveniência das qualidades constantes. Respostas brutas preservadas; não foram descarregados novos ficheiros de sequenciação.

## Decisão e limites

Manter a reanálise v1.22.7 como contagem descritiva dos pares disponibilizados, com esta ressalva explícita. Não remover seis bases automaticamente, não descartar R2 nem interpretar índices de amostra como UMIs sem suporte nos dados. Não promover inferência confirmatória de orientação/protocolo com base apenas no equilíbrio forward/reverse.

Uma análise de sensibilidade limitada comparando R1, R2 e pares pode testar se ambos os mates fornecem alinhamento génico consistente; ainda não foi executada nesta etapa. Não resolveria identidade física nem qual foi o kit. Para fechar proveniência experimental seria necessária documentação original da biblioteca/desmultiplexagem, ficheiros originais ou esclarecimento do depositante. Não foi enviado contacto a ninguém.

As conclusões anteriores de concordância de rótulos e presença de contagens mantêm-se como observações do material analisado. Não atribuir ausência génica a silêncio absoluto. Não há nova validação de acessibilidade ou segurança dos candidatos safe harbor.

## Artefactos

protocol_structure_audit.json, nove XML SRX, nove respostas submitted.tsv, protocol_v1230.py e finish_protocol_v1230.py. Sem alterações aos FASTQs, às contagens, aos filtros ou ao ranking. A recontagem concluída permanece v1.22.7; esta versão é uma auditoria, não uma nova execução do pipeline.
'''
(O/'RELATORIO_PROTOCOLO_2026-09-19.md').write_text(note,encoding='utf-8');(O/'CURRENT.md').write_text(note,encoding='utf-8')
for f in ['protocol_v1230.py','finish_protocol_v1230.py']:shutil.copy2(f,O/f)
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.23.0 — auditoria do protocolo\n\nRecontagem v1.22.7 concluída; contexto v1.22.9 revisto. Kit declarado e estrutura dos FASTQs ainda não reconciliados. Nenhuma alteração de trimming ou ranking.\n\n[Auditoria de protocolo](06_v1.23.0/RELATORIO_PROTOCOLO_2026-09-19.md) · [Contexto RNA](06_v1.22.9/RELATORIO_CONTEXTO_40_GENES_2026-09-19.md).\n',encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.23.0 auditoria protocolo - 2026-09-19';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!warning] v1.23.0 — protocolo ainda não reconciliado\n> [[{name}]]. Auditoria de 1,8 milhões de reads; sem alteração de trimming nem nova recontagem.\n',1),encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.23.0_LIBRARY_PROTOCOL_AUDIT','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print('Protocol audit saved; experimental provenance remains unresolved.')
