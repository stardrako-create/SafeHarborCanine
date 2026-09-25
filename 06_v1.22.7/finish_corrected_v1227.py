from pathlib import Path
import json,csv,hashlib,shutil,datetime
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.7'
review=json.loads((O/'independent_review.json').read_text());assert review['samples']==9 and review['genes']==42309
changes=review['comparison'];qc=list(csv.DictReader((O/'review_sample_QC.tsv').open(),delimiter='\t'));context=list(csv.DictReader((O/'review_focus_gene_context.tsv').open(),delimiter='\t'))
manifest=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'));detail=[];allchanged=set();focus_changed=set()
for r in manifest:
 def load(root):return {v[0]:list(map(int,v[1:])) for l in (root/'samples'/r['SRR']/'STAR_ReadsPerGene.out.tab').read_text().splitlines()[4:] if (v:=l.split('\t'))}
 old=load(P/'06_v1.22.6');new=load(O)
 for g in old:
  if old[g][0]!=new[g][0]:
   allchanged.add(g);detail.append(dict(sample=r['column'],gene=g,old_unstranded=old[g][0],new_unstranded=new[g][0]))
 for row in context:
  if row['sample']==r['column'] and old['gene-'+row['gene']][0]!=new['gene-'+row['gene']][0]:focus_changed.add(row['gene'])
with (O/'existing_gene_count_changes.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['sample','gene','old_unstranded','new_unstranded'],delimiter='\t');w.writeheader();w.writerows(detail)
total=sum(x['added_gene_counts_unstranded'] for x in changes);abschange=sum(x['absolute_count_change_unstranded'] for x in changes)
table='| Amostra | Pares de entrada | Mapeamento único | Pares nos 60 genes recuperados | Genes anteriores alterados |\n|---|---:|---:|---:|---:|\n'
for q,c in zip(qc,changes):
 assert q['sample']==c['sample'];table+=f"| {q['sample']} | {int(q['input_pairs']):,} | {q['unique_pct']} | {c['added_gene_counts_unstranded']:,} | {c['changed_existing_genes_unstranded']} |\n"
focus='| Gene | Intervalo de contagens brutas nas 9 amostras | Intervalo CPM* |\n|---|---:|---:|\n'
for g in dict.fromkeys(r['gene'] for r in context):
 rows=[r for r in context if r['gene']==g];n=[int(r['unstranded_raw_count']) for r in rows];c=[float(r['CPM_gene_assigned_unstranded']) for r in rows];focus+=f'| {g} | {min(n)}–{max(n)} | {min(c):.4f}–{max(c):.4f} |\n'
note=f'''# Recontagem corrigida v1.22.7 — conclusão e revisão, 19 setembro 2026

**Execução computacional concluída e consistência dos outputs verificada nas nove amostras. A validação biológica de safe harbor continua por fazer.** Este estado substitui os checkpoints anteriores de execução.

## Correção e validação

Corrigido exclusivamente o campo source do GTF que continha espaços (`Curated Genomic`), incompatível com a tokenização do STAR. Coordenadas e identificadores mantidos. O novo índice contém os 42.309 genes esperados, incluindo os 60 antes omitidos (61 exões). Referência ROS_Cfam_1.0 e STAR 2.7.11b; esta não é reprodução exata do pipeline original CanFam3.1/HTSeq.

Recontadas as nove amostras, sequencialmente, reutilizando os pares tratados já disponíveis. O worker verificou hashes de referência, índice e inputs tratados, finalização STAR e coerência entre reads pós-fastp, pares de entrada e soma das categorias GeneCounts. Uma revisão separada verificou novamente hashes dos outputs e matrizes, 42.309 IDs por matriz e igualdade célula a célula com cada ficheiro individual nas três orientações. Os certificados MD5 dos downloads tinham sido auditados na v1.22.6; não se alegou novo MD5 integral dos brutos nesta revisão.

## Impacto medido

Os 60 genes recuperados recebem **{total:,} pares atribuídos no total das nove amostras** (contagens não orientadas). Isto não representa novos reads sequenciados: são atribuições possibilitadas pela anotação corrigida. Entre os genes já presentes, mudaram {len(allchanged)} genes distintos, em {len(detail)} combinações gene/amostra; soma de diferenças absolutas = {abschange} contagens. As categorias não atribuídas podem também variar.

{table}

Os dez genes de contexto revistos têm contagens brutas não orientadas {'inalteradas em todas as amostras' if not focus_changed else 'alteradas para: '+', '.join(sorted(focus_changed))}. A correção não autoriza promover os candidatos a safe harbors nem altera diretamente os resultados ATAC. O denominador CPM mudou, pelo que contagem bruta igual não implica CPM exatamente igual.

## Contexto génico corrigido

{focus}

*CPM = contagem não orientada × 1 milhão / soma de contagens génicas não orientadas da mesma amostra. Não assumir igualdade com o denominador da tabela publicada. Os intervalos combinam condições e cães; não são testes de diferenças, equivalência ou estabilidade.

CD247 e CD8B podem agora ser avaliados na anotação ROS e apresentam contagens substanciais; a ausência na tabela CPM publicada não significava ausência de expressão. ANO2, NTF3 e NPNT têm contagens baixas neste conjunto. Isso não demonstra silêncio absoluto, acessibilidade da janela candidata ou segurança de integração. Genes vizinhos expressos, como TBCK e TSEN15, continuam relevantes para avaliar perturbação regulatória, mas expressão por si não demonstra interação com o intervalo.

## Limitações e próximos passos científicos

- Resolver ou manter explicitamente os quatro conflitos nos campos GEO de identidade/condição. Rótulos de biblioteca, título e desenho publicado não são uma confirmação física da identidade. Não efetuar trocas silenciosas.
- A contagem forward/reverse aproximadamente equilibrada é consistente com biblioteca não orientada; reconciliar com o kit SMART-Seq v4 3’ DE declarado e a história de processamento. Não introduzir cortes de reads só pelo nome do kit.
- Q30=100% é pouco informativo neste material com scores constantes de 30 na amostra inspecionada. Não apresentar como qualidade perfeita.
- Comparar expressão com o ficheiro publicado mediante correspondência inequívoca de genes e normalização documentada; não tratar discrepâncias ROS/CanFam3.1 ou STAR/HTSeq como erro automaticamente.
- Qualquer inferência deve respeitar três cães pareados em três condições. Não executar testes de expressão diferencial sobre CPM como se fossem contagens brutas; esta revisão é descritiva.
- RNA génico de CAR-T oferece contexto celular relevante, mas não substitui acessibilidade do locus, estabilidade de expressão do transgene, impacto em genes vizinhos e segurança genómica/funcional. Nenhum locus está validado como safe harbor.

## Espaço e proveniência

Espaço livre no fecho: D {shutil.disk_usage(O).free/1e9:.1f} GB; C {shutil.disk_usage('C:/').free/1e9:.1f} GB. Sem novos downloads; não foi necessário apagar inputs. Ficheiros originais e outputs anteriores preservados para reprodução. O worker terminou e não deve ser reiniciado para repetir esta execução. O acompanhamento desta recontagem será pausado após guardar este relatório.

Evidência: execution_complete.json; index_validated.json; samples/*/validated.json; independent_review.json; comparison_with_v1226.json; existing_gene_count_changes.tsv; review_sample_QC.tsv; review_focus_gene_context.tsv; gene_counts_unstranded.tsv, forward e reverse. Scripts corrected_rna_v1227.py, review_corrected_v1227.py e finish_corrected_v1227.py.
'''
report=O/'RELATORIO_RECONTAGEM_CORRIGIDA_2026-09-19.md';report.write_text(note,encoding='utf-8')
summary='# Estado atual: v1.22.7 concluída e verificada\n\nNove amostras recontadas com anotação corrigida; 42.309 genes. Execução validada; validação biológica de safe harbor pendente.\n\n[Relatório completo](06_v1.22.7/RELATORIO_RECONTAGEM_CORRIGIDA_2026-09-19.md).\n'
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT_before_completion.md');(P/'CURRENT.md').write_text(summary,encoding='utf-8');shutil.copy2(O/'CURRENT.md',O/'running_checkpoint_before_completion.md');(O/'CURRENT.md').write_text(note,encoding='utf-8');shutil.copy2('finish_corrected_v1227.py',O/'finish_corrected_v1227.py')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.7 recontagem concluida - 2026-09-19';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC_before_completion.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!success] v1.22.7 — recontagem corrigida concluída e verificada em 19 setembro\n> [[{name}]]. Nove amostras, 42.309 genes. Substitui estados anteriores de execução; não constitui validação biológica de safe harbor.\n',1),encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.7_CORRECTED_RECOUNT_COMPLETED','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_completion_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(total_recovered_gene_counts=total,distinct_existing_genes_changed=len(allchanged),absolute_existing_change=abschange,focus_genes_changed=sorted(focus_changed),report=str(report)),indent=2))
