# v1.22.7 — recontagem corrigida iniciada em 18 setembro 2026

Novo índice STAR em construção com GTF corrigido. O worker exige os 42.309 genes, hashes da referência e índice, lock exclusivo e validação de outputs. Processa sequencialmente as nove amostras reutilizando trimmed FASTQs existentes, sem downloads nem nova duplicação dos FASTQs. Preserva resultados v1.22.6 para comparação. Cada amostra requer finalização STAR, reconciliação fastp/pares, todas as categorias de contagem e IDs corretos. Matrizes e comparação serão produzidas só após nove amostras validadas.

Espaço inicial: D 160,9 GB; C 61,0 GB. Pausa se D<35 GB ou C<15 GB. Nenhum input foi apagado. Não há novas contagens neste checkpoint. Identidade, orientação e interpretação científica permanecem pendentes.

Estado: 06_v1.22.7/status.json. Script: corrected_rna_v1227.py. Auditoria anterior: 06_v1.22.6/audit_2026-09-18/REVISAO_PIPELINE.md.
