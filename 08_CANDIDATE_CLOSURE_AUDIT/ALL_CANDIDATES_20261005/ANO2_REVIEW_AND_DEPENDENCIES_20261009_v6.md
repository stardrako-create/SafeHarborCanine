# ANO2-39708 — revisão consolidada de09/10/2026

**Candidato adicional com suporte aos filtros implementados de locus; aprovação integral da auditoria/intervenção ainda não emitida. Nenhuma validação biológica.**

ROS_Cfam_1.0: **NC_051831.1:39708074–39709074**, BED0,end-exclusive,1000bp. Projeção UU:NC_049248.1:6741835–6742835, orientação reversa. É outro locus relativamente a w01/w11; pertence à região da reserva w03 e não deve contar como um locus independente de w03.

| Verificação | Resultado | Limite |
|---|---|---|
| Origem e ATAC | Região SHIP convergente original62427bp; percentis regional72,70 e20kb65,97 ≥55 | Sangue público, não célulasCAR-T |
| RNA e distâncias ROS | Sem overlap consultado; mRNA5′ mais próximo102152bp,risco880858bp,miRNA1273355bp | Anotações, não isolamento funcional |
| RNA UU | Zero overlap de genes/pseudogenes,exões,CDS/RNAs; mRNA e RNA5′ mais próximos:NTF3,101951bp | NCBI AR106; fronteiras anotadas, não TSS medidos |
| Risco/miRNA UU | GAPDH879877bp;miRNA1270499bp;zero hits associados dentro de300kb | 2547/2654 símbolos relacionados após aliases/GeneIDs;107não resolvidos, nunca tratados como ausência de risco |
| Repetidos e mapeamento | RepeatMasker:0bp na janela; ROS1000M,NM0,MAPQ60,sem alternativas/XA/SA | Biblioteca/source e referência pesquisados |
| Conservação |1000basesROS; máximo média50bp1,55246<6,5, após retirada de espécies duplicadas | Suporte desigual2–194espéciesACGT/base; reanálise do mesmo alinhamento |
| Regulação/cromatina | EpiC:1000/1000bases anotadas em11fontes,zero estados promotor/enhancer;TAD sem hits consultados | Proxies públicos; não isolamento regulatório emCAR-T |
| SV/alelos locais | Nenhum hit no catálogo141810registos/150BND;5indelsPASS normalizados;63120mers de máscara intactos em32combinações | Não exclusão deSV individuais nem32haplótipos observados; fase não inventada |
| Metilação | RRBS150/1000bp,média48,73% | Não demonstra hipometilação de toda a janela |

## Especificidade de referência: duas hipóteses, não guias aprovados

| Hipótese | Protospacer ROS BED0 | Sequência5′→3′ | PAM |
|---|---|---|---|
| ANO2_TTTV23_1 |39708151–39708174,+ |TATGAGGATTTGTATCAAACACA |TTTC |
| ANO2_TTTV23_2 |39708436–39708459,+ |CCTCGGCTGTGTAATTCCTGGAG |TTTC |

As duas têm zero outros hits **no modelo definido**, tanto em ROS completo como UU completo:23nt;PAMsTTTV/CTTV/TCTV/TTCV;até3substituições sem gap,ou exatamente um bulge de1base com até2substituições,ambas as orientações. O índice UU foi conferido contra todos os2197identificadores/comprimentos da FASTA,total2481983352bases; hashes dos seis ficheiros e do genoma guardados. Cada alvo pretendido foi encontrado uma vez na orientação UU prevista.

PAM/protospacer locais permanecem intactos nas32combinações sintéticas dos cinco indelsPASS. A sensibilidade a uma única mudança identificou19e4sítios hipotéticos, respetivamente:23projeções UU contíguas foram confrontadas com os catálogos;nenhum SNP recuperado criou hits no modelo. Os três indels recuperados e todas as40combinações sintéticas dos alelos nesses sítios também não criaram hits.

**Limites:** PAMs/variantes de enzima não cobertos,múltiplos bulges,alelos fora dos23sítios,genomas individuais e lacunas/bases ambíguas não estão excluídos. Bulges sobrepostos ao alvo pretendido foram separados; não foi demonstrada ausência de cortes deslocados locais. Comparação de referências e reanálise das mesmas leituras/alinhamentos não são confirmação biológica independente. O painelNGG anterior permanece preservado, com correspondências semelhantes; não foi criado um novo veto por ter hits a3mismatches nem escolhida uma nuclease para obter uma aprovação favorável.

## Stop concreto da aprovação integral

A **nuclease e variante realmente previstas não foram indicadas**. É necessário esse parâmetro para avaliar compatibilidade doPAM/crRNA, aplicar o modelo de especificidade da enzima real e fechar a seleção do guia. Os plasmídeosDrosophila históricos não definem o sistema canino. O cão experimental não foi exigido para completar estes scans.

Os107símbolos de risco sem correspondência permanecem uma limitação de cobertura de anotação, explicitamente documentada. Não se declara ausência de todos os riscos ou prova funcional de segurança. Genótipo do futuro cão,cassete e ensaios de integração,expressão,genes vizinhos e funçãoCAR-T são dependências individuais/experimentais separadas. Não são resultados já obtidos.

## Expansão independente e fecho da etapa executada

Conservação concluída após revisão taxonómica:1059=0,50206;1266=2,21532;1900=0,46692,1000bases cada,abaixo de6,5. Os três gaps estão fora da geometria SHIP original;1059 também tem ATAC do gap insuficiente e projeções alélicas pendentes;1059/1900 não são convergentes. Estes resultados não foram usados para reabilitar vetos ou mudar gates. Os inputs,contagens e outputs anteriores continuam preservados.

HAL,scoring,pesquisa UU e revisões encadeadas terminaram. A assertividade inicial do índice falhou porque o resumo contém descriçõesFASTA completas e oFAI usa o primeiro identificador; a versão2 normalizou apenas essa representação e confirmou igualdade de todos os identificadores/comprimentos. Tentativas e falhas anteriores preservadas. Não há worker ativo desta etapa.

## Recibos e proveniência

No diretório `../INDEPENDENT_LOCI_20261008/`:
- `ANO2_two_guides_full_UU_20261009/results.json`,com alinhamentos brutos e comandos;
- `full_UU_index_integrity_20261009_v2/receipt.json` e `full_UU_sequence_provenance_20261009.json`;
- `ANO2_TTTV23_offsite_UU_and_catalogue_20261009/results.json`,VCFs brutos e projeção;
- `ANO2_offsite_three_indels_20261009/receipt.json` e `ANO2_offsite_catalogue_combinations_20261009/receipt.json`;
- `ANO2_UU_annotation_crosscheck_20261009/receipt.json` e `ANO2_UU_risk_alias_GeneID_crosscheck_20261009.json`;
- `remote_conservation_taxa_review_20261009/receipt.json`;
- os recibos originais de ANO2 e versões anteriores do relatório,referenciados na versão5.

Anotação:[NCBI AR106](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Canis_lupus_familiaris/106/),GGF com MD5 oficial e SHA256 conferidos. CURRENT e Obsidian atualizados. MemPalace permanece indisponível (SSE404),com fila local preservada; não alegar sincronização. Nenhum contacto externo ou publicação.
