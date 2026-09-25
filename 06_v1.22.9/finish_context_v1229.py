from pathlib import Path
import csv,json,shutil
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.9';s=json.loads((O/'context_summary.json').read_text());assert s['mapped']==40 and not s['unresolved']
rows=list(csv.DictReader((O/'all40_overview.tsv').open(),delimiter='\t'));table='| Gene | Ausente na tabela publicada | Amostras com contagem >0 | Contagens brutas min–max | CPM min–max |\n|---|---|---:|---:|---:|\n'
for r in rows:table+=f"| {r['gene']} | {'sim' if r['absent_published_CPM']=='True' else 'não'} | {r['nonzero_samples']}/9 | {r['min_raw']}–{r['max_raw']} | {float(r['min_CPM']):.4f}–{float(r['max_CPM']):.4f} |\n"
note=f'''# v1.22.9 — contexto CAR-T dos 40 genes previamente selecionados

19 setembro 2026. **Os 40 genes têm correspondência exata e única na anotação corrigida; os 15 ausentes da matriz publicada podem agora ser avaliados nesta reanálise.** Não são 15 genes recuperados pelo bug dos 60 genes: são problemas diferentes. A primeira ausência refere-se à tabela publicada; a segunda era a ingestão do GTF pela v1.22.6.

## Âmbito e validação

Mantida a seleção definida na v1.22.4: dez genes distintos mais próximos de cada uma das três âncoras declaradas, mais marcadores/genes predefinidos. Não se selecionaram genes pela expressão observada. Mantidas as âncoras históricas, sem as apresentar como coordenadas finais de integração. Os nomes históricos dos loci não substituem a lista efetiva de vizinhos e distâncias em 06_v1.22.4/local_gene_context.tsv.

Usada a matriz não orientada corrigida v1.22.7, cuja consistência foi verificada. Correspondência pelo gene_name único do índice ROS, sem aliases forçados. Produzidos 360 registos gene/amostra, 120 resumos gene/condição e 360 comparações descritivas dentro do mesmo rótulo de cão. Os totais foram verificados e os hashes de inputs guardados. A v1.22.8 já demonstrou melhor correspondência 9/9 com os rótulos publicados; isto não confirma identidade física.

CPM é contagem bruta × 1 milhão / soma dos pares atribuídos a genes na mesma biblioteca. É normalização descritiva por tamanho, sem correção de composição; não são valores TPM, nem necessariamente o CPM do autor. Comparações usam três cães reportados, B/E/M, cada um com controlo, CAR simples e CAR duplo.

## Resultados

{table}

CD247 e CD8B apresentam contagens em todas as amostras apesar de ausentes da tabela CPM publicada. Isso impede interpretar aquela ausência como ausência de expressão. ANO2 (0–4 pares) e NTF3 (0–4) têm suporte de contagem muito baixo; ratios relativos neste domínio são frágeis. LOC119876429 apresenta 0–32 pares, e LOC119872513 zero em todas as amostras. KCNA5 e KCNA6 também têm zero pares atribuídos. Contagem zero é uma observação sob esta biblioteca, anotação e regra de atribuição — não prova de silêncio absoluto.

## Comparações pareadas

within_donor_descriptive_ratios.tsv contém CAR simples/controlo, CAR duplo/controlo e duplo/simples por cão. O log2-ratio só existe quando ambos os CPM são estritamente positivos. Não se adicionou pseudocount para criar razões artificiais nas contagens zero. Valores ausentes por zero e por anotação têm estados explícitos. Os ratios de genes com pouquíssimos pares devem ser lidos com os respetivos numeradores/denominadores, não como efeitos fiáveis.

condition_descriptive_summary.tsv apresenta mediana e amplitude entre três cães em cada condição. Não foram calculados p-values, declarada estabilidade/equivalência, nem usados estes resumos como exclusão automática. Ausência de significância num eventual teste futuro também não provaria neutralidade.

## Implicações para os candidatos

Esta análise fecha a lacuna de contexto RNA dos 15 símbolos ausentes da tabela publicada. Não fecha o problema ATAC nem demonstra expressão estável de um transgene. Genes vizinhos com contagens substanciais continuam a exigir avaliação de possível perturbação regulatória; baixa expressão de outro vizinho não torna o intervalo automaticamente seguro.

Não foi alterado o ranking, nenhum filtro de candidatos foi afrouxado e nenhuma âncora foi promovida a alvo final. Continuam relevantes os conflitos regulatórios, cobertura e identidade celular do ATAC, incertezas de variantes populacionais e a falta de validação funcional pós-integração documentados nas auditorias anteriores.

## Pendências reais

Reconciliação do kit/protocolo e dos quatro conflitos source/tissue; correspondência génica mais ampla se for necessária inferência transcriptómica; normalização apropriada e modelo pareado se houver uma pergunta de expressão diferencial; evidência celular e funcional de acessibilidade, expressão persistente e segurança da integração. Nenhum contacto externo foi enviado nem nova execução pesada iniciada nesta etapa.

Artefactos: gene_mapping.tsv; all40_gene_sample_context.tsv; all40_overview.tsv; condition_descriptive_summary.tsv; within_donor_descriptive_ratios.tsv; context_summary.json. Scripts context_v1229.py e finish_context_v1229.py. A recontagem concluída continua a ser a v1.22.7.
'''
(O/'RELATORIO_CONTEXTO_40_GENES_2026-09-19.md').write_text(note,encoding='utf-8');(O/'CURRENT.md').write_text(note,encoding='utf-8')
for n in ['context_v1229.py','finish_context_v1229.py']:shutil.copy2(n,O/n)
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.22.9 — contexto RNA dos 40 genes revisto\n\n40/40 genes com correspondência única, incluindo 15 ausentes na matriz publicada. Recontagem v1.22.7 validada; concordância v1.22.8 revista. Não constitui validação biológica de safe harbor.\n\n[Contexto detalhado](06_v1.22.9/RELATORIO_CONTEXTO_40_GENES_2026-09-19.md) · [Concordância](06_v1.22.8/RELATORIO_CONCORDANCIA_2026-09-19.md).\n',encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.9 contexto 40 genes - 2026-09-19';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!info] v1.22.9 — contexto dos 40 genes\n> [[{name}]]. 15 símbolos ausentes da tabela publicada agora avaliáveis; RNA não valida safe harbor.\n',1),encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.9_ALL40_CAR_T_RNA_CONTEXT','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print('40-gene context documented in project, Obsidian and MemPalace.')
