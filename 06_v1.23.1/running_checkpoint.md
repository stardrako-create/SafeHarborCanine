# v1.23.1 — teste limitado R1/R2/pares em execução

19 setembro 2026. Iniciados 27 alinhamentos sequenciais: nove amostras × paired, R1 e R2. Para cada amostra são usados os mesmos primeiros 100.000 pares dos FASTQs tratados da v1.22.6, com o índice corrigido v1.22.7. Não foram adicionados filtros, alterados dados completos ou descarregadas novas sequências. Esta comparação controla a seleção dos reads entre os três modos.

Objetivo: verificar contribuição de ambos os mates para mapeamento/atribuição génica, como análise de sensibilidade perante a discrepância entre nome do kit e estrutura dos FASTQs. Não identifica o kit nem resolve a identidade física dos cães.

Proteções: lock exclusivo; verificação de nomes dos pares, tamanho sequência/qualidade, hashes dos subsets, ID dos 42.309 genes, contagem de 100.000 unidades de entrada por modo, reconciliação das categorias GeneCounts e marcador de finalização STAR. Os recibos guardam parâmetros e hashes. O índice corrigido é reutilizado; a tabela geneInfo é verificada contra o recibo anterior. Não se recalcularam neste teste os hashes dos grandes binários do índice.

Inputs single-end são reads e paired são pares: as taxas devem ser interpretadas com esse denominador explícito. A correlação de contagens usa genes não zero em pelo menos um dos dois modos. O subset inicial não é amostragem aleatória e é insuficiente para interpretar zeros de genes raros ou estabilidade entre condições.

Uso adicional estimado abaixo de 1 GB em subsets e outputs; sem BAM. Reservas mínimas de 20 GB em D e 15 GB em C, verificadas por amostra. Não é necessário apagar os FASTQs nesta etapa. Nenhum resultado novo concluído neste checkpoint.

Estado status.json; saída final execution_complete.json, mate_alignment_QC.tsv e mate_count_concordance.tsv. O acompanhamento automático foi ativado e deverá pausar após a revisão final. Última recontagem completa validada permanece v1.22.7.
