# Controlos: evidência integrada

As 20 janelas foram revistas. Esta tabela junta as medições e mantém as lacunas visíveis; não aprova um painel final. Variantes parciais só descrevem trechos exatos. Conservação conta apenas posições com score e bases comparativas de outra espécie.

| Controlo | Repetidos bp/1000 | 20mers únicos/981 | Regulatório bp/1000 | Tecidos com alerta | Variantes: janela inteira / parcial | Conservação com suporte bp/1000 |
|---|---:|---:|---:|---|---|---:|
| bg1k_2054 | 389 | 696 | 1000 | LI, ST | — / 28 | 1000 |
| bg1k_0567 | 223 | 630 | 980 | nenhum nas bases avaliadas | 45 / — | 980 |
| bg1k_2848 | 0 | 787 | 1000 | nenhum nas bases avaliadas | 24 / — | 1000 |
| bg1k_0045 | 365 | 528 | 879 | LI | — / 67 | 879 |
| bg1k_2865 | 305 | 928 | 946 | nenhum nas bases avaliadas | 16 / — | 946 |
| bg1k_0442 | 128 | 918 | 1000 | nenhum nas bases avaliadas | 35 / — | 1000 |
| bg1k_2901 | 802 | 914 | 1000 | nenhum nas bases avaliadas | — / 10 | 1000 |
| bg1k_0202 | 57 | 705 | 1000 | nenhum nas bases avaliadas | — / 50 | 1000 |
| bg1k_0334 | 310 | 959 | 991 | nenhum nas bases avaliadas | — / 22 | 991 |
| bg1k_2545 | 21 | 952 | 979 | LU | — / 31 | 979 |
| bg1k_0228 | 0 | 717 | 1000 | nenhum nas bases avaliadas | — / 21 | 1000 |
| bg1k_1958 | 491 | 51 | 1000 | nenhum nas bases avaliadas | — / 21 | 1000 |
| bg1k_2686 | 492 | 316 | 1000 | nenhum nas bases avaliadas | — / 37 | 1000 |
| bg1k_1023 | 195 | 446 | 987 | nenhum nas bases avaliadas | 47 / — | 987 |
| bg1k_0023 | 742 | 63 | 81 | nenhum nas bases avaliadas | não resolvido | 62 |
| bg1k_0878 | 335 | 550 | 1000 | nenhum nas bases avaliadas | 93 / — | 807 |
| bg1k_1650 | 225 | 169 | 935 | MG | — / 47 | 935 |
| bg1k_2228 | 160 | 922 | 1000 | PA | 22 / — | 1000 |
| bg1k_2867 | 81 | 812 | 590 | nenhum nas bases avaliadas | não resolvido | 590 |
| bg1k_0117 | 1000 | 15 | 1000 | nenhum nas bases avaliadas | 37 / — | 0 |

## Restrições por janela

- bg1k_2054: TAD_proxy_risk_gene; EpiC_promoter_enhancer_in_observed_bases; UU_whole_window_not_exact
- bg1k_0567: TAD_proxy_risk_gene; TAD_proxy_missing_assignment; TAD_boundary_overlap; regulatory_projection_incomplete; conservation_supported_coverage_incomplete
- bg1k_2848: Sem estes alertas nas camadas avaliadas; seleção final e validação experimental pendentes.
- bg1k_0045: EpiC_promoter_enhancer_in_observed_bases; regulatory_projection_incomplete; UU_whole_window_not_exact; conservation_supported_coverage_incomplete
- bg1k_2865: regulatory_projection_incomplete; UU_multiple_reported_alignments; conservation_supported_coverage_incomplete
- bg1k_0442: Sem estes alertas nas camadas avaliadas; seleção final e validação experimental pendentes.
- bg1k_2901: UU_whole_window_not_exact
- bg1k_0202: UU_whole_window_not_exact
- bg1k_0334: regulatory_projection_incomplete; UU_whole_window_not_exact; conservation_supported_coverage_incomplete
- bg1k_2545: TAD_proxy_risk_gene; EpiC_promoter_enhancer_in_observed_bases; regulatory_projection_incomplete; UU_whole_window_not_exact; conservation_supported_coverage_incomplete
- bg1k_0228: UU_whole_window_not_exact
- bg1k_1958: TAD_proxy_risk_gene; UU_whole_window_not_exact
- bg1k_2686: UU_whole_window_not_exact
- bg1k_1023: regulatory_projection_incomplete; conservation_supported_coverage_incomplete
- bg1k_0023: regulatory_projection_incomplete; UU_multiple_reported_alignments; UU_whole_window_not_exact; conservation_supported_coverage_incomplete
- bg1k_0878: conservation_supported_coverage_incomplete
- bg1k_1650: TAD_proxy_missing_assignment; EpiC_promoter_enhancer_in_observed_bases; regulatory_projection_incomplete; UU_whole_window_not_exact; conservation_supported_coverage_incomplete
- bg1k_2228: EpiC_promoter_enhancer_in_observed_bases
- bg1k_2867: regulatory_projection_incomplete; UU_multiple_reported_alignments; UU_whole_window_not_exact; conservation_supported_coverage_incomplete
- bg1k_0117: UU_multiple_reported_alignments; conservation_supported_coverage_incomplete

## O que continua por resolver

Bases sem projeção e alelos que atravessam diferenças/gaps; regiões com múltiplos alinhamentos; origem faseada dos marcadores de deleção. Qualquer avaliação por janela inteira deve manter estas lacunas. w01 continua prioritário e w11 mantém a dúvida estrutural já documentada. A escolha da nuclease, cassete e cães dadores continua necessária para o desenho específico.

Os ficheiros de origem e hashes estão em review.json. Outputs anteriores foram preservados.