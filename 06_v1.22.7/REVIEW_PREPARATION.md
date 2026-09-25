# v1.22.7 — revisão preparada, execução em curso

Verificada uma única instância do worker e STAR a construir o índice corrigido. Sem falha registada. D aproximadamente 161 GB e C 61 GB livres no checkpoint.

Preparado review_corrected_v1227.py: só produz resultados após execution_complete.json, verifica hashes de matrizes e outputs individuais, reconcilia todas as contagens e exporta QC e contexto de dez genes com CPM explicitamente definido como contagem por milhão de pares atribuídos a genes. Não confundir esse denominador com o CPM publicado. O script foi executado antes da conclusão e recusou corretamente produzir resultados.

O GEO SOFT local GSE247355_family.soft.gz declara SMART-Seq v4 3’ DE Kit. A formulação publicada sobre filtragem Phred >30 é ambígua e não foi traduzida automaticamente numa opção de processamento. As qualidades constantes de 30 no FASTQ e o equilíbrio forward/reverse não resolvem a história de preparação/processamento da biblioteca.

Guardado read_structure_check.json: primeiros 10.000 reads de cada um dos 18 FASTQs; métricas de hexâmeros iniciais, GC e frequência de pelo menos 16 T nas posições 7–26. Esta amostragem inicial é uma verificação exploratória, não aleatória nem confirmação do protocolo. Nenhum trimming adicional foi aplicado com base no nome do kit.

O acompanhamento automático já está ativo. Após a recontagem, executar review_corrected_v1227.py e interpretar comparison_with_v1226.json antes de promover resultados. Identidade amostral, adequação inferencial e validação de safe harbor continuam em aberto.
