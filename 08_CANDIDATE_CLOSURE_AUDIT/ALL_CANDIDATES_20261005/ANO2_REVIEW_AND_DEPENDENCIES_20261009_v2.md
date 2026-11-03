# ANO2-39708: revisão computacional e decisões pendentes

ROS_Cfam_1.0, BED0: **NC_051831.1:39708074–39709074** (1000 bp). Coordenadas inclusivas a partir de1:39708075–39709074.

É uma alternativa adicional a w01 e w11. Pertence à região ANO2/NTF3 onde está a reserva w03; estas janelas não são dois loci independentes.

**Resultado de locus:** os filtros implementados de seleção de locus têm suporte positivo nos dados interrogados, com os mesmos limiares anteriores. **Resultado da auditoria completa:** ainda sem aprovação da intervenção/guia. **Resultado biológico:** não validado. Não publicar nem comunicar este documento como demonstração de segurança de integração ou de expressão CAR-T.

| Critério de locus | Resultado revisto | Alcance |
|---|---|---|
| Origem SHIP | Região convergente ANO2/NTF3,62427bp | Região original entre as461; não extensão escolhida depois |
| ATAC regional | Percentil72,70 ≥55 | Background3000 com distribuição original de comprimentos |
| ATAC local20kb | Percentil65,97 ≥55 | Track de sangue, não demonstra acessibilidade em célulasCAR-T |
| Repetidos |0bp na janela na nova execução com flancos10kb | Biblioteca dog; flancos com hits, janela semhit |
| Conservação |1000basesROS; máximo da média50bp1,56426 <6,5 | Projeção porbase;11bases exclusivasCanFam3 e4diferenças explícitas; sensibilidadeROS concordante |
| Mapeamento |1000M,NM0,MAPQ60; semXA/SA/alternativos | BWA aoROS completo; regra de locus anterior satisfeita |
| Anotações de RNA/regulação | Sem overlap nas fontes consultadas | GFF/transcritos/TSS anotados,smallRNA/lncRNA e mapas externos; ausência de hit não é isolamento funcional |
| Distâncias anotadas | mRNA5′ mais próximo102152bp; gene de risco880858bp; miRNA1273355bp | Limiares anteriores50/300/300kb; extremidade anotada, não TSS medido |
| EpiC |1000/1000bases anotadas emcada uma de11fontes;0bp promotor/enhancer | Interseção porbase, não envelope1011bp; tecidos públicos |
| TAD | Sem fronteira/hit de risco nas escalas consultadas | Proxy deHi-C, não contactos medidos emCAR-T |
| SV | Nenhumhit no rastreio141810registos, incluindo150BND | Catálogo consultado; não exclusão deSV individuais |
| Projeção/alelos |UU reverso1000M,NM1, discrepância reconciliada;5indelsPASS normalizados | REF/haplótipos conferidos;3símbolos* tratados com união conservadora e coocorrência genotípica, sem assumir fasecis |
| Robustez de sequência |63120mers fora máscaraREF intactos em32combinações dos5indels;102protospacers+PAMNGG intactos | Combinações sintéticas, não32haplótipos observados; não genótipos do cão experimental |

## Especificidade: resultados, sem aprovação de guias

A pesquisa exata tem918/98120mers únicos;631 ficam fora da máscaraREF aplicada. Todos901100mers são únicos,44fora da máscara. O scan de semelhanças contraROS completo revalidou1135513alinhamentos para2127queries de3loci.

O painel exploratório20nt+NGG contém102sequências emANO2.41não têm outra correspondência comNGG/PAMindefinido até2mismatches, mas todas têm correspondências a3mismatches. Na pesquisa adicional de umbulge e até2substituições,9dessas41não têm hits fora do locus até1substituição; todas têm hits a2substituições. Os9perfis foram anotados com genes: têm correspondências em regiões génicas. São candidatos a avaliação, não guias aprovados nem evidência de edição nesses sítios.

A definição NGG é uma hipótese técnica de comparação. Não foi escolhida a enzima do laboratório. Uma correspondência com3mismatches não foi transformada num novo veto ao locus; também não foi descartada como biologicamente irrelevante. PAMalternativos, variante exata da nuclease, múltiplosbulges, genomasalternativos e atividade real não estão fechados.

O painel alternativo23nt+TTTV terminou para80sequências em3loci. EmANO2,5/14não têm hits no modelo canónico; ao alargar aPAMs com uma substituiçãoT>C,2/14continuam semhits no modelo definido. É um diagnóstico para apoiar a escolha. [A documentação do fabricante](https://www.neb.com/en-us/products/m0653-engen-lba-cas12a-cpf1) descreveTTTV; [o comprimento de targeting](https://www.neb.com/faqs/how-do-i-design-a-guide-rna-for-use-with-engen-lba-cas12a) permite18–24nt.23nt é um parâmetro exploratório explícito. [Um estudo primário sobreLbCas12a](https://academic.oup.com/nar/article/48/7/3722/5763098) mostra que PAMnão canónicos também podem importar; um scan apenascanónico nunca demonstra ausência de todos os off-targets da enzima real.

## O que falta e o que não bloqueia a análise de referência

1. Especificar nuclease e variante/PAM realmente previstos; depois fechar a avaliação dos guias e dos sítios de integração adequados à cassete. A pergunta está pendente. O histórico dos plasmídeosDrosophila não define o sistema canino.
2. Não assumir um haplótipo individual nem uma fasecis inexistente. A união conservadora resolve robustez local das sequências testadas, sem inventar confirmação alélica.
3. Genótipo do futuro cão e ensaios de integração, expressão persistente, efeitos emgenesvizinhos e funçãoCAR-T são dependências individuais/experimentais. Não são necessários para continuar a auditoria computacional sobreROS; não se apresentam como já validados.

RRBS cobre150/1000bp, média48,73%, máximo12cães porlocal. Não há evidência de hipometilação da janela inteira. Zero inícios nas célulasT rasas disponíveis é insuficiência de suporte direto, não veto demonstrado. Nenhum limiar novo foi criado para promover ou excluir esta janela.

## Recibos e reprodução

Os resultados e scripts novos estão em `../INDEPENDENT_LOCI_20261008/`. Ver `ANO2_closure_checkpoint_20261009.json`, `whole_window_ROS_mappability_20261009/receipt.json`, `ANO2_EpiC_basewise_review_20261009.json`, `ANO2_original_regional_ATAC_20261009.json`, `ANO2_conservation_identity_20261009_v3/receipt.json`, `ANO2_ROS_model_sensitivity_20261009/receipt.json`, `ANO2_indel_union_robustness_20261009.json`, `ANO2_hypothetical_NGG_PAM_indel_robustness_20261009.json`, `hypothetical_NGG_single_bulge_interpretation_20261009.json` e `ANO2_hypothetical_NGG_gene_context_20261009.json`. Hashes/fontes e falhas anteriores preservados. Não alterar contagens históricas para concordarem com um novo estado.

MemPalace: MCPindisponível(SSE404); atualização pendente na fila local. Obsidian eCURRENT atualizados. Nenhuma publicação ou contacto externo foi efetuado nesta revisão.


## Adendo: duas hipóteses de sequência com resultado computacional favorável

| Hipótese | Protospacer ROS BED0 | Strand | PAM | Resultado do modelo |
|---|---|---|---|---|
| ANO2_TTTV23_1 |NC_051831.1:39708151–39708174 |+ |TTTC |0hits fora do locus no modelo definido |
| ANO2_TTTV23_2 |NC_051831.1:39708436–39708459 |+ |TTTC |0hits fora do locus no modelo definido |

Modelo explícito:23nt; PAMsTTTV/CTTV/TCTV/TTCV; alinhamentosungapped até3substituições ou exatamenteumbulge de1base até2substituições; referênciaROS pesquisável. Foram testadas as duas orientações. Protospacer ePAM mantêm-se intactos nas32combinações sintéticas dos5indelsPASS locais. A construção e a pesquisa foram conferidas com exemplos de22/23/24bases contra enumeração exaustiva. Recibos e sequências em `ANO2_TTTV23_REFERENCE_MODEL_CANDIDATES_20261009.tsv/json` no diretório de análise.

**Isto não define a nuclease real nem aprova guias.** O modelo não abrange todas as especificidades dePAM de todas as variantes da enzima, múltiplosbulges, variantes emoutrosloci ou genomasindividuais. Basesambíguas/lacunas doassembly são informação emfalta, não prova deausência deoff-targets. Os alinhamentos comumbulge que sobrepõem o intervalo pretendido foram separados do rastreio de outros loci; não se demonstrou ausência decortes deslocados localmente. Zero neste modelo não significa zerooff-targets biológicos.

## Retoma após interrupção

09/10/2026: o utilizador informou que teve dedesligar. A tentativaHALv2 ficou comestado`running`,mas semprocessos presentes; a causa precisa no sistemaoperativo não foi estabelecida. OMAFparcial foi preservado semvalidação. Iniciado umúnico workerHALv3 destacado (launcherWindows27748,wrapperLinux292), usando cacheexistente e novosoutputs. Aumentado apenasocap operacional de transferência nova para512MiB; pedidosHTTPS mantêm validação206/Content-Range/ETag e gravação durável. Score só decorre após recibodeextração`completed`. Não alterados limiares científicos. Estado em `verified_HAL_v3_chain_20261009/status.json` e `remote_HAL_verified_20261009_v3/status.json`.
