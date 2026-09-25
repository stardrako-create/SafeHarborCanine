# Auditoria do pipeline v1.22.6 — 18 setembro 2026

As nove amostras terminaram, mas a reanálise fica **provisória: foi identificado um defeito na anotação ingerida pelo STAR**. A execução concluída não equivale a validação científica.

## Verificações realizadas

Os nove logs STAR contêm ALL DONE e as consolas finished successfully. As contagens de entrada concordam exatamente com metade dos reads pós-fastp; as quatro categorias especiais e os genes somam o total de pares em cada orientação. As três matrizes agregadas concordam célula a célula com os nove outputs. São 42.249 linhas génicas.

Os 18 FASTQs têm tamanho concordante com ENA e certificados prévios com MD5 esperado=observado. Não se repetiu nesta auditoria o hash integral dos 63,9 GB; a auditoria não deve ser descrita como nova verificação integral dos FASTQs. Foram inspecionados 10.000 reads por ficheiro.

Mapeamento único: 89.13–93.54%. Atribuição génica não orientada: 74.87–79.85%. A fração forward entre contagens orientadas é 50,04–50,21%, consistente com biblioteca não orientada; confirmar protocolo antes de análise inferencial. Não escolher uma coluna orientada apenas por ser ligeiramente maior.

## Achado prioritário: perda silenciosa de anotação

O campo source `Curated Genomic` contém espaço. O STAR 2.7.11b lê os primeiros campos por whitespace e deixa de reconhecer essas linhas como exon. O nosso validador anterior usava tabs, pelo que não detetou a incompatibilidade. Foram ignorados 61 exões, afetando 60 genes: 60 completamente ausentes e 0 parcialmente afetados. A diferença de 60 genes entre GTF e índice é exatamente explicada pelos genes que perderam todos os exões.

Preparado `ROS_release106.STAR_compatible.gtf`, alterando apenas espaços do campo source para underscores; coordenadas, atributos e IDs mantidos. Validação com a tokenização do STAR recupera 42309 genes e 1062532 exões. O índice original e as contagens não foram sobrescritos. É necessário reconstruir o índice e recontar as nove amostras antes de promover os resultados a definitivos. Mesmo genes sem exões diretamente perdidos podem ter atribuição alterada em sobreposições; não basta acrescentar 60 linhas de zeros.

Fonte técnica: [parser GTF do STAR 2.7.11b](https://github.com/alexdobin/STAR/blob/2.7.11b/source/GTF.cpp), leitura dos campos por stream e reconhecimento da feature exon.

## Qualidade de sequenciação e interpretação

Todos os 18 FASTQs amostrados mostram apenas Phred 30 (27 milhões de bases no total), e o fastp reporta Q30=100% em todas as amostras. Isto limita o valor informativo do QC de qualidade de bases; não demonstra sequenciação perfeita. Scores simplificados são uma possibilidade, não proveniência confirmada destes ficheiros. [Documentação NCBI sobre qualidade simplificada](https://www.ncbi.nlm.nih.gov/sra/docs/sra-data-formats/).

fastp removeu pequenas quantidades de reads mesmo com filtros genéricos desativados (incluindo adaptadores/dímeros e reads curtos); o texto anterior não deve sugerir retenção absoluta. Os totais pós-fastp foram conciliados com STAR. As contagens orientadas estão preservadas conforme o [manual STAR](https://github.com/alexdobin/STAR/blob/2.7.11b/doc/STARmanual.pdf).

CD247 e CD8B, ausentes da tabela CPM publicada, têm contagens abundantes nesta anotação. ANO2/NTF3 têm contagens muito baixas. Isto esclarece ausência na tabela versus ausência de expressão, mas os valores permanecem provisórios até à correção do índice. RNA génico não prova acessibilidade nem neutralidade dos intervalos safe harbor.

## Robustez operacional e espaço

Há aproximadamente 161.1 GB livres em D: neste checkpoint. Não foram apagados downloads: são úteis para corrigir esta execução e o espaço atual permite preservá-los. O pipeline processa alinhamentos sequencialmente, mas descarrega até seis ficheiros e conserva raw+trimmed; não implementa ainda a estratégia ponta a ponta de uma amostra de cada vez com libertação após validação solicitada pelo utilizador.

Problemas de recuperação: a fila aceita presença de ficheiros como evidência de conclusão sem assinatura de inputs/parâmetros; reinícios reinicializam o histórico de comandos; não há lock exclusivo; um antigo queue_failure.json permanece ao lado do marcador de sucesso. Estes problemas não invalidam por si os nove outputs aqui conferidos, mas devem ser corrigidos antes de reutilizar a fila. Esta auditoria não alterou retroativamente logs nem apagou o erro histórico.

Antes de qualquer limpeza: conclusão STAR comprovada, integridade das contagens, concordância com referência e parâmetros, relatórios QC e checksums duráveis. A unidade operacional é o par R1/R2 da mesma amostra. A limpeza deve abranger apenas os FASTQs dessa amostra e ser registada; conservar ficheiros necessários à análise corrigida. O orçamento deve contar D: e também C: onde reside o WSL/tmp.

## Pendências científicas

Recontagem com anotação compatível; confirmar orientação/protocolo; avaliar concordância amostral com CPM publicado e identidade perante os quatro conflitos GEO; respeitar desenho de três cães pareados; avaliar transformação/normalização e sensibilidade antes de inferência. Não há ainda demonstração de safe harbor em CAR-T caninas: continuam em falta evidência de cromatina no contexto relevante e validação funcional/segurança de integração.

Evidência reprodutível: audit.json, sample_QC.tsv, STAR_annotation_parser_impact.json e script audit_run_v1226.py. Nenhum novo ranking de candidatos foi aprovado nesta auditoria.
