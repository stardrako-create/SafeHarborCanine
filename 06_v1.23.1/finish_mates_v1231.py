from pathlib import Path
import json,shutil
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.23.1';s=json.loads((O/'independent_review.json').read_text());assert s['status']=='independent_review_passed' and len(s['validated_alignments'])==27
table='| Modo | Mapeamento único (%) | Atribuição génica (%) |\n|---|---:|---:|\n'
for mode,v in s['mode_ranges'].items():table+=f"| {mode} | {v['unique_pct'][0]:.2f}–{v['unique_pct'][1]:.2f} | {v['assigned_pct'][0]:.3f}–{v['assigned_pct'][1]:.3f} |\n"
ct='\n'.join(f'- {k}: {v[0]:.4f}–{v[1]:.4f}.' for k,v in s['correlation_ranges'].items())
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file())
note=f'''# v1.23.1 — teste dos mates concluído e revisto

19 setembro 2026. **27/27 alinhamentos concluídos e verificados independentemente. R1 e R2 fornecem alinhamento genómico e atribuição génica substanciais nas nove amostras.** O teste não sustenta descartar R2 como se fosse composto apenas por índices/poly(T) nos ficheiros tratados analisados.

## Desenho e validação

Por amostra, os mesmos primeiros 100.000 pares tratados foram alinhados em três modos: paired, apenas R1 e apenas R2. Nove amostras; 900.000 pares distintos entre as amostras, reutilizados nos três modos. Não são 2,7 milhões de pares independentes. Índice corrigido v1.22.7, STAR 2.7.11b. Não foram alterados parâmetros do pipeline completo nem contagens v1.22.7.

O revisor verificou os 27 hashes de contagens, hashes dos subsets contra recibos, finalização STAR, 42.309 IDs esperados, soma das categorias GeneCounts=100.000 unidades em cada orientação, concordância das tabelas com os logs e recalculou as 27 correlações. Todas as verificações passaram. A unidade single-end é um read; no modo paired é um par. As taxas abaixo usam esses denominadores; não somar contagens R1+R2 como moléculas independentes.

## Resultados nas nove amostras

{table}

Correlações de Spearman de contagens não orientadas, restringindo cada comparação a genes não zero em pelo menos um dos modos:

{ct}

Os pares têm maior taxa de mapeamento único e atribuição génica que cada mate isolado neste teste. A proximidade de R1 e R2 reforça que ambos contêm sequência utilizável. Não há fundamento neste resultado para remover R2 nem cortar seis bases automaticamente. A concordância é parcialmente esperada porque os mates vêm dos mesmos fragmentos e não constitui replicação biológica independente.

## O que fica por resolver

O teste fornece suporte empírico ao uso operacional dos pares disponibilizados, mas não identifica o kit original, não reconstrói desmultiplexagem/processamento anterior e não resolve os quatro conflitos de metadados. A discrepância com a estrutura descrita no manual continua documentada na v1.23.0. O presente teste usou reads tratados; a auditoria anterior avaliou a estrutura dos brutos disponíveis.

Amostragem inicial determinística, não aleatória, e profundidade de apenas 100.000 pares por amostra. Não interpretar zeros de genes raros, diferenças pequenas entre modos ou resultados dos loci candidatos como testes de estabilidade, equivalência ou segurança. Não foram gerados p-values. A contagem completa v1.22.7 e o contexto dos 40 genes v1.22.9 permanecem as referências para descrição de expressão.

Não foi demonstrada acessibilidade em CAR-T dos intervalos candidatos, persistência do transgene, ausência de perturbação regulatória nem segurança da integração. Nenhum candidato foi promovido a safe harbor.

## Espaço e fecho

Pasta desta análise: {size/1e9:.3f} GB, abaixo da estimativa de 1 GB. Espaço livre observado no fecho: D {shutil.disk_usage(O).free/1e9:.1f} GB e C {shutil.disk_usage('C:/').free/1e9:.1f} GB. A queda total de espaço em D durante o intervalo excede o tamanho desta pasta; a origem do restante consumo não foi determinada e não foi atribuída a este teste. Não foram apagados ficheiros de outros processos nem inputs não validados. O worker terminou.

O acompanhamento será pausado após este fecho. Artefactos: independent_review.json, execution_complete.json, mate_alignment_QC.tsv, mate_count_concordance.tsv e recibos por amostra/modo. Scripts mate_sensitivity_v1231.py, review_mates_v1231.py e finish_mates_v1231.py.
'''
(O/'RELATORIO_MATES_2026-09-19.md').write_text(note,encoding='utf-8');shutil.copy2(O/'CURRENT.md',O/'running_checkpoint.md');(O/'CURRENT.md').write_text(note,encoding='utf-8')
for f in ['review_mates_v1231.py','finish_mates_v1231.py']:shutil.copy2(f,O/f)
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT_before_completion.md');(P/'CURRENT.md').write_text('# Estado atual: v1.23.1 concluída e revista\n\n27/27 alinhamentos do teste R1/R2/pares verificados. Ambos os mates têm suporte genómico e génico substancial. Proveniência experimental ainda não reconciliada; não constitui validação biológica de safe harbor.\n\n[Relatório dos mates](06_v1.23.1/RELATORIO_MATES_2026-09-19.md). Recontagem completa de referência: v1.22.7.\n',encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.23.1 teste mates concluido - 2026-09-19';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC_before_completion.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!info] v1.23.1 — teste dos mates concluído\n> [[{name}]]. 27/27 verificados; R1 e R2 informativos. Não confirma o kit nem valida safe harbor.\n',1),encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.23.1_MATE_SENSITIVITY_COMPLETED','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_completion_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print('Mate sensitivity report saved to project, Obsidian and MemPalace.')
