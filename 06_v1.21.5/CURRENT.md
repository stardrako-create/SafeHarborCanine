# Checkpoint v1.21.5 — limpeza técnica

2026-09-14. Relatório de apresentação adiado por instrução do utilizador. Este é apenas um checkpoint de execução.

Concluído: leitura direta de 76 cães em 3 janelas locais e 3 regiões (456 medições), sem erro de leitura nesses intervalos; comparação 71/76 com os mesmos pesos QC e sem confundir diferenças de ganho; bootstrap emparelhado de 2 000 reamostragens dos 71 cães em janelas fixas; sensibilidade a pesos/QC/gate; reconstrução do consenso de 76 cães (39 votos, merge 75 bp) e scoring completo coerente. Os 461 scores, estados e vetos ATAC são iguais aos da v1.21.3.

Isso preserva três passes das verificações legadas, não três loci livres de todas as evidências agora disponíveis. A anotação independente EpiC Dog acrescenta sinais regulatórios em cada uma das três regiões. Nas janelas previamente escolhidas de 1 kb: ANO2 tem 800 bp na união de estados de enhancer em pulmão/glândula mamária; NPNT tem 662 bp em rim; LOC119876429 tem 0 bp de promotor/enhancer nos 11 tecidos. As janelas de 1 kb com alguma sobreposição são 78/247, 122/227 e 39/225, respetivamente. Nenhum destes tecidos equivale a células T purificadas; ausência de anotação não demonstra neutralidade.

Fontes EpiC: https://doi.org/10.1126/sciadv.ade3399 e https://github.com/snu-cdrc/dog-reference-epigenome . Os BEDs foram fixados por commit, com integridade conferida por tamanho e SHA-1 de blob Git; provenance em external_sources/epic_sources.json. Classes amplas promotor 1–4 e enhancer 5–7 conferidas no script original Figure 6C_1_re-category_chromatin_states.sh. Não importámos a classificação humana.

ATAC local: 64/76 cães acima da própria média genómica em ANO2 e NPNT, 68/76 em LOC119876429. Frequência de pico por cão: 0, 9 e 3, respetivamente. As cinco amostras adicionais reduzem a média gated-CPM local em 4,83%, 5,04% e 6,14%. Bootstrap das diferenças entre candidatos inclui zero em todos os pares. Isto não corrige a seleção prévia dos máximos nem demonstra estabilidade em CAR-T.

Remapeamento local: 308/309 fragmentos de 100/150/250 bp passam o proxy definido (MAPQ≥30, alinhamento exato esperado, sem alternativa reportada pelo BWA). A exceção é um fragmento de 100 bp em LOC119876429. Não é ensaio exaustivo de mappability nem análise off-target de gRNA.

Dog10K SV: as três janelas alinham integralmente e sem alternativa reportada em UU_Cfam_GSD_1.0 (MAPQ60; NM 6/1/0 em ANO2/LOC/NPNT). ANO2 e LOC não têm registos SV sobrepostos no catálogo consultado. NPNT sobrepõe uma deleção agregada de ~28,9 Mb, AF≈0,1064%, quatro genótipos heterozigóticos, três com FT PASS. A variante é um flag a rever, não uma deleção demonstrada nestes cães nem um veto automático. Fonte: https://kiddlabshare.med.umich.edu/dog10K/Manta-SV_2022-03-28/SV-genotype-v2.merge.agg_only.08032022.vcf.gz . SNPs/indels ainda não avaliados.

QC: só 2/76 têm FRiP≥0,2 no cálculo local contra picos individuais. Não foi declarado QC absoluto aprovado. O TSS enrichment local é uma média em janelas largas e não a métrica de pico TSS ENCODE. A sensibilidade aos 33 cães acima de ambas as medianas QC altera valores relativos e troca LOC/NPNT; não estabelece um vencedor robusto.

Falta: avaliação específica de T cells/CAR-T; QC e contexto de expressão adequados; verificação de SNPs/indels e da chamada SV grande; modelo de conservação independente (o existente é piloto); tratamento coerente da evidência EpiC na decisão local; reconstrução ATAC genome-wide com recuperação completa. Leitura local recuperada não equivale a rebuild global. Não foi escolhido um ponto de inserção, removido um veto de segurança ou declarado um safe harbor validado.

Dados operacionais: cohort_and_support.tsv, per_dog_local_atac.tsv, fixed_window_bootstrap.tsv, paired_bootstrap_differences.tsv, technical_sensitivity.tsv, local_remap_tiles.tsv, dog10k_lookup.json, NPNT_SV_review.json, epic_local_states.tsv, local_windows_regulatory_rechecked.tsv, regulatory_recheck_summary.json, validation.json. Código e comandos preservados na pasta.
