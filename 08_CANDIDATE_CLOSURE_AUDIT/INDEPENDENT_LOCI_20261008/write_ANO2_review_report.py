from pathlib import Path
import json,hashlib,datetime
N=Path(__file__).resolve().parent;A=N.parent/'ALL_CANDIDATES_20261005';out=A/'ANO2_REVIEW_AND_DEPENDENCIES_20261009.md';assert not out.exists()
text='''# ANO2-39708: revisão computacional e decisões pendentes

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

Um painel alternativo23nt+TTTV está em execução, como diagnóstico para apoiar a escolha. [A documentação do fabricante](https://www.neb.com/en-us/products/m0653-engen-lba-cas12a-cpf1) descreveTTTV; [o comprimento de targeting](https://www.neb.com/faqs/how-do-i-design-a-guide-rna-for-use-with-engen-lba-cas12a) permite18–24nt.23nt é um parâmetro exploratório explícito. [Um estudo primário sobreLbCas12a](https://academic.oup.com/nar/article/48/7/3722/5763098) mostra que PAMnão canónicos também podem importar; um scan apenascanónico nunca demonstra ausência de todos os off-targets da enzima real.

## O que falta e o que não bloqueia a análise de referência

1. Especificar nuclease e variante/PAM realmente previstos; depois fechar a avaliação dos guias e dos sítios de integração adequados à cassete. A pergunta está pendente. O histórico dos plasmídeosDrosophila não define o sistema canino.
2. Não assumir um haplótipo individual nem uma fasecis inexistente. A união conservadora resolve robustez local das sequências testadas, sem inventar confirmação alélica.
3. Genótipo do futuro cão e ensaios de integração, expressão persistente, efeitos emgenesvizinhos e funçãoCAR-T são dependências individuais/experimentais. Não são necessários para continuar a auditoria computacional sobreROS; não se apresentam como já validados.

RRBS cobre150/1000bp, média48,73%, máximo12cães porlocal. Não há evidência de hipometilação da janela inteira. Zero inícios nas célulasT rasas disponíveis é insuficiência de suporte direto, não veto demonstrado. Nenhum limiar novo foi criado para promover ou excluir esta janela.

## Recibos e reprodução

Os resultados e scripts novos estão em `../INDEPENDENT_LOCI_20261008/`. Ver `ANO2_closure_checkpoint_20261009.json`, `whole_window_ROS_mappability_20261009/receipt.json`, `ANO2_EpiC_basewise_review_20261009.json`, `ANO2_original_regional_ATAC_20261009.json`, `ANO2_conservation_identity_20261009_v3/receipt.json`, `ANO2_ROS_model_sensitivity_20261009/receipt.json`, `ANO2_indel_union_robustness_20261009.json`, `ANO2_hypothetical_NGG_PAM_indel_robustness_20261009.json`, `hypothetical_NGG_single_bulge_interpretation_20261009.json` e `ANO2_hypothetical_NGG_gene_context_20261009.json`. Hashes/fontes e falhas anteriores preservados. Não alterar contagens históricas para concordarem com um novo estado.

MemPalace: MCPindisponível(SSE404); atualização pendente na fila local. Obsidian eCURRENT atualizados. Nenhuma publicação ou contacto externo foi efetuado nesta revisão.
'''
out.write_text(text,encoding='utf-8')
manifestfiles=[N/'ANO2_closure_checkpoint_20261009.json',N/'ANO2_hypothetical_NGG_PAM_indel_robustness_20261009.json',N/'hypothetical_NGG_single_bulge_interpretation_20261009.json',N/'ANO2_hypothetical_NGG_gene_context_20261009.json',out]
manifest=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in manifestfiles],complete_audit_approval=False,biological_validation=False)
(N/'ANO2_review_delivery_manifest_20261009.json').write_text(json.dumps(manifest,indent=2));print(str(out))
