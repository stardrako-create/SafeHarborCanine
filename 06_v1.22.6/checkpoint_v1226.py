from pathlib import Path
import json,shutil,datetime
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.6';a=json.loads((O/'reference/annotation_validation.json').read_text());assert a['status']=='passed';d=json.loads((O/'ENA_verified_download_manifest.json').read_text());assert len(d['files'])==18
genes=set(a['gene_ids']);required=['ANO2','NTF3','CD247','CD8B'];coverage={g:sorted(x for x in genes if x in [g,'gene-'+g]) for g in required};assert all(coverage.values())
(O/'critical_gene_annotation_presence.json').write_text(json.dumps(coverage,indent=2))
note='''# Safe Harbor CAR-T — v1.22.6: reprocessamento bruto em execução

Iniciado em 16 de setembro de 2026. **Este é um checkpoint de execução, não um resultado concluído.** O último contexto RNA concluído é a v1.22.4, com auditoria de fontes v1.22.5. Nenhuma contagem nova ou conclusão biológica é anunciada aqui.

## Trabalho iniciado

- Manifesto ENA atualizado para os nove runs: 18 FASTQs paired-end, tamanhos e MD5 oficiais, estratégia RNA-Seq e ligações SRX/SRR verificadas.
- Downloads retomáveis iniciados, com até seis ficheiros em paralelo. A transferência está lenta e poderá exigir muitas horas. Um ficheiro só se torna elegível para análise depois de tamanho e MD5 concordarem com o manifesto. Ficheiros parciais não entram na contagem.
- Ambiente isolado `cart_rna` instalado: STAR 2.7.11b e gffread 0.12.7. Exportação explícita dos pacotes guardada. O fastp disponível é 1.3.6 e as opções usadas foram verificadas no executável.
- Referência local ROS_Cfam_1.0 GCF_014441545.1 e GFF NCBI Release106 fixados por SHA256. Não é a reprodução exata da análise original CanFam3.1.
- Conversão GTF validada: 42.309 identificadores génicos com exões, 1.062.532 features de exões. Coordenadas contra o FASTA e associação inequívoca de cada transcrito a um gene verificadas.
- A conversão gffread omitia gene_id em features de transcritos ligados diretamente a genes/pseudogenes. Foram completadas 7.910 features exclusivamente quando transcript_id coincide com o ID original de gene/pseudogene e o nome concorda; as coordenadas não foram alteradas. A conversão bruta e a reparação estão preservadas.
- ANO2, NTF3, CD247 e CD8B têm identificador na anotação de contagem preparada. Isto confirma que podem ser avaliados nesta reanálise; não significa que já tenham reads ou expressão detetada.
- Construção do índice STAR iniciada, com seis threads, SA esparso e limite de geração de 28 GB, tendo em conta os 31 GiB de RAM disponíveis no WSL.

## Fila de processamento

Foi iniciado um processo de fila que aguarda índice concluído e pares de FASTQ verificados. Processa uma amostra de cada vez. O fastp fará deteção/remoção de adaptadores paired-end e relatórios QC, mantendo desativados os filtros genéricos de qualidade/comprimento e a remoção automática de poly-G. Não há deduplicação UMI nem alegação de contagem molecular. O impacto e o protocolo de biblioteca terão de ser revistos nos relatórios QC.

O STAR fará alinhamento genómico anotado, single-pass, e GeneCounts. Não serão produzidos BAMs nesta execução para limitar armazenamento; reads brutos, reads tratados, logs, QC e contagens são conservados. Isto limita a revisão posterior de alinhamentos e identidade, que poderá exigir uma execução dedicada.

Serão preservadas as três colunas de contagem — não orientada, forward e reverse — sem escolher silenciosamente a orientação. Só após as nove execuções serão agregadas matrizes. Mesmo nessa altura, o estado será execução concluída com revisão científica pendente: qualidade de mapeamento, strandedness, identidade amostral, correspondência de genes e adequação a análise diferencial têm de ser avaliadas.

A fila aborta perante erro de processamento ou falha registada de download; não contorna validação. Tem limite operacional de sete dias. O sistema de downloads tenta retomar transferências interrompidas. Os processos dependem de o computador/WSL continuar disponível; não são um serviço remoto nem uma garantia de conclusão.

## Como acompanhar

`download_status.json`, logs `*.download.log`, `queue_status.json`, `reference/index_console.log`, `reference/index_execution.json` quando existir, e eventual `queue_failure.json`. A dimensão de um ficheiro ainda aberto pode não refletir imediatamente o progresso mostrado pelo curl; a confirmação final usa tamanho e MD5 após fecho.

Saídas finais esperadas: `gene_counts_unstranded.tsv`, `gene_counts_forward.tsv`, `gene_counts_reverse.tsv` e `raw_reanalysis_execution_complete.json`. A existência deste último indica execução, não validação biológica ou aprovação de safe harbor.

Scripts: `download_cart_v1226.py`, `prepare_star_v1226.py`, `run_cart_queue_v1226.py`. Os processos foram iniciados nesta sessão; antes de retomar noutro turno, verificar processos e estados para evitar duplicar downloads, índices ou alinhamentos. A anotação e preparação estão verificadas; downloads, índice e alinhamentos ainda não concluídos neste checkpoint.
'''
(O/'CURRENT.md').write_text(note,encoding='utf-8')
for name in ['download_cart_v1226.py','prepare_star_v1226.py','run_cart_queue_v1226.py','checkpoint_v1226.py']:shutil.copy2(Path(name),O/name)
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.6 RNA bruto em execucao - 2026-09-16';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!important] Em execução: v1.22.6 — RNA bruto\n> [[{name}]]. Downloads com MD5, índice STAR e fila iniciados. Ainda sem novas contagens. Último relatório de resultados concluídos: v1.22.5.\n',1),encoding='utf-8')
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.22.6 em execução\n\n[Checkpoint da execução](06_v1.22.6/CURRENT.md) · [Último relatório concluído](06_v1.22.5/RELATORIO_ATUALIZADO_2026-09-16.md).\n\nDownloads, índice e fila de reprocessamento iniciados. Não anunciar contagens ou conclusão antes de verificar os estados e outputs.\n',encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.6_RAW_RNA_RUNNING','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_start_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print('Running checkpoint saved; no completion claimed.');print(coverage)
