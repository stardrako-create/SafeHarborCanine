exec(open('audit_run_v1226.py',encoding='utf-8').read().split('checks=[]; quality=[]')[0])
from mem_local import call
data=json.loads((A/'audit.json').read_text()); missing=set(data['missing_annotation_genes']); bad=collections.Counter(); genes=collections.Counter(); total=collections.Counter()
source=O/'reference/ROS_release106.gtf'; fixed=A/'ROS_release106.STAR_compatible.gtf'; fixes=0
with source.open() as inp,fixed.open('w') as out:
 for line in inp:
  v=line.rstrip('\n').split('\t')
  if len(v)==9:
   if v[2]=='exon':
    gid=re.search(r'gene_id "([^"]+)"',v[8])[1];total[gid]+=1
    if len(v[1].split())!=1:bad[v[1]]+=1;genes[gid]+=1
   if len(v[1].split())!=1:v[1]='_'.join(v[1].split());line='\t'.join(v)+'\n';fixes+=1
  out.write(line)
fully={g for g,n in genes.items() if n==total[g]}; assert fully==missing
# Validate the repaired file with the same whitespace tokenization STAR uses.
parsed=set(); nexon=0
for line in fixed.open():
 v=line.split()
 if len(v)>2 and v[2]=='exon': parsed.add(re.search(r'gene_id "([^"]+)"',line)[1]);nexon+=1
assert parsed==set(json.loads((O/'reference/annotation_validation.json').read_text())['gene_ids'])
assert nexon==json.loads((O/'reference/annotation_validation.json').read_text())['exon_features']
impact=dict(ignored_exons_by_source=dict(bad),affected_gene_exons=dict(genes),fully_missing_genes=sorted(fully),partially_affected_genes=sorted(set(genes)-fully),changed_source_fields=fixes,repaired_exons=nexon,repaired_genes=len(parsed),repaired_sha256=hashlib.sha256(fixed.read_bytes()).hexdigest(),existing_index_unchanged=True)
(A/'STAR_annotation_parser_impact.json').write_text(json.dumps(impact,indent=2))
note=f'''# Auditoria do pipeline v1.22.6 — 18 setembro 2026

As nove amostras terminaram, mas a reanálise fica **provisória: foi identificado um defeito na anotação ingerida pelo STAR**. A execução concluída não equivale a validação científica.

## Verificações realizadas

Os nove logs STAR contêm ALL DONE e as consolas finished successfully. As contagens de entrada concordam exatamente com metade dos reads pós-fastp; as quatro categorias especiais e os genes somam o total de pares em cada orientação. As três matrizes agregadas concordam célula a célula com os nove outputs. São 42.249 linhas génicas.

Os 18 FASTQs têm tamanho concordante com ENA e certificados prévios com MD5 esperado=observado. Não se repetiu nesta auditoria o hash integral dos 63,9 GB; a auditoria não deve ser descrita como nova verificação integral dos FASTQs. Foram inspecionados 10.000 reads por ficheiro.

Mapeamento único: {min(r['unique_pct'] for r in data['samples']):.2f}–{max(r['unique_pct'] for r in data['samples']):.2f}%. Atribuição génica não orientada: {min(r['assigned_pct'] for r in data['samples']):.2f}–{max(r['assigned_pct'] for r in data['samples']):.2f}%. A fração forward entre contagens orientadas é 50,04–50,21%, consistente com biblioteca não orientada; confirmar protocolo antes de análise inferencial. Não escolher uma coluna orientada apenas por ser ligeiramente maior.

## Achado prioritário: perda silenciosa de anotação

O campo source `Curated Genomic` contém espaço. O STAR 2.7.11b lê os primeiros campos por whitespace e deixa de reconhecer essas linhas como exon. O nosso validador anterior usava tabs, pelo que não detetou a incompatibilidade. Foram ignorados {sum(bad.values())} exões, afetando {len(genes)} genes: {len(fully)} completamente ausentes e {len(set(genes)-fully)} parcialmente afetados. A diferença de 60 genes entre GTF e índice é exatamente explicada pelos genes que perderam todos os exões.

Preparado `ROS_release106.STAR_compatible.gtf`, alterando apenas espaços do campo source para underscores; coordenadas, atributos e IDs mantidos. Validação com a tokenização do STAR recupera {len(parsed)} genes e {nexon} exões. O índice original e as contagens não foram sobrescritos. É necessário reconstruir o índice e recontar as nove amostras antes de promover os resultados a definitivos. Mesmo genes sem exões diretamente perdidos podem ter atribuição alterada em sobreposições; não basta acrescentar 60 linhas de zeros.

Fonte técnica: [parser GTF do STAR 2.7.11b](https://github.com/alexdobin/STAR/blob/2.7.11b/source/GTF.cpp), leitura dos campos por stream e reconhecimento da feature exon.

## Qualidade de sequenciação e interpretação

Todos os 18 FASTQs amostrados mostram apenas Phred 30 (27 milhões de bases no total), e o fastp reporta Q30=100% em todas as amostras. Isto limita o valor informativo do QC de qualidade de bases; não demonstra sequenciação perfeita. Scores simplificados são uma possibilidade, não proveniência confirmada destes ficheiros. [Documentação NCBI sobre qualidade simplificada](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/).

fastp removeu pequenas quantidades de reads mesmo com filtros genéricos desativados (incluindo adaptadores/dímeros e reads curtos); o texto anterior não deve sugerir retenção absoluta. Os totais pós-fastp foram conciliados com STAR. As contagens orientadas estão preservadas conforme o [manual STAR](https://github.com/alexdobin/STAR/blob/2.7.11b/doc/STARmanual.pdf).

CD247 e CD8B, ausentes da tabela CPM publicada, têm contagens abundantes nesta anotação. ANO2/NTF3 têm contagens muito baixas. Isto esclarece ausência na tabela versus ausência de expressão, mas os valores permanecem provisórios até à correção do índice. RNA génico não prova acessibilidade nem neutralidade dos intervalos safe harbor.

## Robustez operacional e espaço

Há aproximadamente {data['disk_free_bytes']/1e9:.1f} GB livres em D: neste checkpoint. Não foram apagados downloads: são úteis para corrigir esta execução e o espaço atual permite preservá-los. O pipeline processa alinhamentos sequencialmente, mas descarrega até seis ficheiros e conserva raw+trimmed; não implementa ainda a estratégia ponta a ponta de uma amostra de cada vez com libertação após validação solicitada pelo utilizador.

Problemas de recuperação: a fila aceita presença de ficheiros como evidência de conclusão sem assinatura de inputs/parâmetros; reinícios reinicializam o histórico de comandos; não há lock exclusivo; um antigo queue_failure.json permanece ao lado do marcador de sucesso. Estes problemas não invalidam por si os nove outputs aqui conferidos, mas devem ser corrigidos antes de reutilizar a fila. Esta auditoria não alterou retroativamente logs nem apagou o erro histórico.

Antes de qualquer limpeza: conclusão STAR comprovada, integridade das contagens, concordância com referência e parâmetros, relatórios QC e checksums duráveis. A unidade operacional é o par R1/R2 da mesma amostra. A limpeza deve abranger apenas os FASTQs dessa amostra e ser registada; conservar ficheiros necessários à análise corrigida. O orçamento deve contar D: e também C: onde reside o WSL/tmp.

## Pendências científicas

Recontagem com anotação compatível; confirmar orientação/protocolo; avaliar concordância amostral com CPM publicado e identidade perante os quatro conflitos GEO; respeitar desenho de três cães pareados; avaliar transformação/normalização e sensibilidade antes de inferência. Não há ainda demonstração de safe harbor em CAR-T caninas: continuam em falta evidência de cromatina no contexto relevante e validação funcional/segurança de integração.

Evidência reprodutível: audit.json, sample_QC.tsv, STAR_annotation_parser_impact.json e script audit_run_v1226.py. Nenhum novo ranking de candidatos foi aprovado nesta auditoria.
'''
(A/'REVISAO_PIPELINE.md').write_text(note,encoding='utf-8');shutil.copy2('audit_run_v1226.py',A/'audit_run_v1226.py');shutil.copy2('finish_pipeline_audit.py',A/'finish_pipeline_audit.py')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');dest=notes/'Safe Harbor CAR-T - Auditoria pipeline v1.22.6 - 2026-09-18.md';dest.write_text(note,encoding='utf-8')
(P/'CURRENT.md').write_text('# Estado atual: v1.22.6 executada; correção de anotação necessária\n\nNove amostras contadas. Auditoria de 18 setembro identificou exões ignorados pelo parser STAR devido a espaços no campo source. Resultados provisórios.\n\n[Auditoria detalhada](06_v1.22.6/audit_2026-09-18/REVISAO_PIPELINE.md).\n',encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';moc.write_text(m.replace(anchor,anchor+'\n\n> [!warning] Auditoria 18 setembro — nove amostras concluídas; anotação STAR requer correção.\n> [[Safe Harbor CAR-T - Auditoria pipeline v1.22.6 - 2026-09-18]]. Este estado substitui o checkpoint de execução de 16 setembro.\n',1),encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.6_PIPELINE_AUDIT_20260918','content':note,'source_file':str(dest),'added_by':'Codex'}); assert not r.get('error') and not r.get('result',{}).get('isError');(A/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(impact,indent=2))
