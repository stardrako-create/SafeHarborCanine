# 07_integracao — Cross-Layer Integration and Ranking

Combines the per-locus safe-harbor score (`05_SHIP/`) with gene/element
context (`06_GEG-SH/`) into a single ranked candidate list, applying
exclusion rules (e.g. proximity to oncogenes/tumor suppressors, essential
genes, TAD boundaries) and producing the shortlist that feeds
`08_candidatos_finais/`.

Not yet populated — depends on `05_SHIP/` and `06_GEG-SH/`.
