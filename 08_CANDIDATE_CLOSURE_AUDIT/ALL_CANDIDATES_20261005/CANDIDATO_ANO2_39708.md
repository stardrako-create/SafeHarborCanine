# ANO2-39708 — candidato adicional para avaliação

Data: 08/10/2026. Estado: CANDIDATO ADICIONAL PROPOSTO; aprovação computacional completa ainda não demonstrada. Sem validação biológica. Não substituir este estado por PASS.

## Identidade
ROS_Cfam_1.0 NC_051831.1:39708074–39709074 (BED0, 1000bp). Em coordenadas1-based inclusivas:39708075–39709074. Região ANO2/NTF3, distinta das regiões de w01 e w11. É uma janela diferente da reserva w03 dentro da mesma região ANO2/NTF3; não contar ambas como dois loci independentes.
UU: NC_049248.1:6741835–6742835, reverso,1000M,MAPQ60,NM1. CanFam6:chr27:7264044–7265044, duas rotas concordantes na fonte histórica. A identidadeCanFam6 não implica identidadeUU.

## Evidência revista
- Sem sobreposição com genes, transcritos, miRNA, lncRNA/smallRNA ou elementos regulatórios externos nas anotações interrogadas.
- Início de mRNA anotado mais próximo: NTF3,102152bp. Limite anotado, nãoTSS medido.
- Gene mais próximo do catálogo de risco:GAPDH,880858bp; a categoria é do catálogo, não uma nova conclusão sobre risco.
- miRNA mais próximo:1273355bp.
- Sem sobreposição com fronteiraTAD; todas as escalas de intervalos entre fronteiras consultadas atribuídas, sem gene do catálogo de risco nesses intervalos. Proxy deHi-C, não contactos específicos deCAR-T.
- NenhumhitSV na consulta integral141810registos, incluindo150ALT BND, no envelopeUU. Não excluiSV individual.
- Zero sobreposição repetitiva e promotor/enhancerEpiC na evidência preexistente. Sem nova execuçãoRepeatMasker nesta revisão.
- ATACsangue percentil90,37; evidência histórica do gate min/mediana/max21/35/43cães. Não é acessibilidadeCAR-T.
- RRBS15%observado, média48,73%, máximo12cães porlocal. Não demonstra hipometilação emtodaajanela.
- UnicidadeROS:918/981segmentos20bp únicos;631 permanecem semsobreposição com a máscaraREF aplicada. Todos901segmentos100bp únicos;44 semsobreposição com amáscara. Não são guias, não incluem mismatches/bulges/PAM.
- VarianteROS G>A nooffset555:1700/3974alelos,42,778%; mantida namáscara. ContagensreconciliadasAC/AN, semgenótipoindividual.
- Trêssímbolos* coincidem porposição comdeleçõesPASS2/6bp; associação porhaplótipo não demonstrada. Mantidos comoalertas; não contar comoSNPs independentes nemdescartar.

## Motivo para propor esta janela
Oferece região adicional a w01 com ausência dos hitsSV deNPNT/TBCK, semrepetidos anotados nesta janela e menor metilaçãoobservada que a reservaw03(88,27%em10%). Épriorização do conjunto deevidência, não ranking estatístico nemalegação de superioridadefuncional.

## Pendências que impedem aprovação completa
Concluir contexto dosalelos simbólicos/efeitosindel e avaliação deespecificidade dependente danuclease. As camadas regulatórias públicas não demonstram isolamento funcional. Zero iníciosT nos dadosrasos disponíveis não prova locus fechado, mas tambémnão demonstra acessibilidadeT. Integração, expressão persistente, efeitosnosgenesvizinhos e função celular requeremexperiências. Não afrouxar nenhum destes requisitos para cumprir ameta.

## Fontes locais
UU_projection/projection_review.json; all_RNA_context.json; transcript_and_boundary_review.tsv; TAD_multiscale_review.json; direct_annotation_review_20261008.json; SV_all_windows/status.json; masked_specificity_20261008.json; SNP_reference_rebase_counts_20261008.json; star_and_reference_context_20261006.json; UU_projection/inputs.json. Scripts associados e resultados anteriores preservados.

## Adendo 08/10/2026
Nova execução RepeatMasker dog com flancos10kb confirmou0bp anotados nesta janela (INDEPENDENT_LOCI_20261008/repeatmasker_two_20261008). GateATAC20kb do código legado passouP65,9667 (ATAC_narrow_gate_20261008.json). Estes resultados atualizam as frases históricas acima, sem alterar a avaliação biológica. As dependências experimentais e de escolha de nuclease são distintas das lacunas da auditoria computacional de locus; não se deve apresentar uma exigência experimental como teste computacional já executado ou falhado.

## Adendo 09/10/2026 — contexto genotípico dos símbolos*
Nova consulta completa de genótipos em ANO2_star_genotype_context_20261009: cada um dos três alertas tem um único portador, também portador da deleção próxima correspondente. Não existe fase explícita suficiente para afirmar ligação cis. Mantidos alertas e máscara; associação compatível, sem aprovação integral ou confirmação independente. Outputs anteriores intactos.
## Adendo 09/10/2026 — conservação, contexto ATAC e indels
Conservação recalculada a partir do MAF retido:1000basesROS, máximo50bp1,56426<6,5. Projeção única porbase comsalto explícito de11basesCanFam3; quatrodiferenças pontuais todasna máscaraREF. Sensibilidade comDOGreferênciaROS nascolunas existentes mantém1,56426. Não seafirma alinhamentoindependente.
Gate regionaloriginal62427bp recomputado:P72,70>=55; gate local20kbP65,97>=55. CincoindelsPASS não sobrepostos:631segmentos20bp fora máscara ficamintactos emtodas32combinaçõessintéticas. Isto resolve robustez dessas sequências perante aunião destesindels, não afase nemtodosgenomasalternativos. Especificidade até3mismatches emtodoROS emexecução; ainda semaprovaçãointegral/biológica.

## 09/10/2026 — mapeamento e especificidade concluídos; fecho de guia pendente
ANO2-39708: BWA do intervalo1000bp=1000M/NM0/MAPQ60, sem XA/SA/alternativos; regra de mapeamento existente satisfeita. EpiC revisto por1000basesROS únicas em11fontes: cobertura integral emcadafonte e0bp promotor/enhancer; não prova isolamento funcional.
ScanBowtie completo:2127queries/1135513alinhamentos revalidados, todos exatos únicos, nenhum semnearmatches<=3 semPAM. Painel exploratório20nt+NGG: ANO2 tem102sequências,96 semoutras NGG/indefinidasaté1mismatch e41até2;0até3.1266:97/93/46/1;1900:82/73/39/0 (total/semoutras<=1/<=2/<=3). Não são guiasaprovados nemcritério novo de veto ao locus. Nuclease real perguntada ao utilizador, mantendo análisesindependentes emcurso.
Checkpoint ANO2_closure_checkpoint_20261009.json contém hashes/fontes, limites e separação locus/guia/experiências. Os filtros de locus têm suporte revisto; não se declara aprovaçãointegral daauditoria enquanto especificidade da intervenção não estiver fechada. Conservação de1059/1266/1900 continua comworkerúnico e cacheverificada; nenhumduplicate lançado. MemPalace SSE404: atualização apenasnafilalocal.


## 09/10/2026 — retoma após interrupção e painéis de especificidade aprofundados
HALv2 estava semworker e comrecibo`running`desatualizado após desligar;MAFparcial preservado,nãovalidado. HALv3 único iniciado destacado, novosoutputs e cacheexistente; cap512MiB novos, mesmoscritérios. LauncherWindows27748/wrapperLinux292; verificarverified_HAL_v3_chain_20261009/status.json e remote_HAL_verified_20261009_v3/status.json antes derelançar.
ANO2:102protospacers+PAMNGG intactos em32combinaçõesindel;9perfisNGG revistos contraanotaçõesgénicas e nenhumaprovado. Painel23ntTTTV:14hipótesesANO2;5semhitsno modelocanónico;2mantêmzero aoalargar paraTTTV/CTTV/TCTV/TTCV,até3substituições ungapped ou1bulge até2substituições. Estes2alvos/PAM estãointactos em32combinações. Basesambíguas/genomasalternativos, atividadeenucleasereal nãocobertos;nãosão guiasebiologiavalidados. Nucleaserealcontinua perguntada, semresposta nestaexecução.
Relatório novoANO2_REVIEW_AND_DEPENDENCIES_20261009_v2.md preservaov1eassuas fontes. Checkpoints anterioresnãoalterados;nãodeclararaprovaçãointegraldaauditoria.MemPalaceSSE404: filalocal,nãosync.


## 09/10/2026 — controlo de conservação com espécies duplicadas
MAFretidoANO2:4105pares espécie/coluna duplicados,47comACGT discordante. Inverter ordem altera65scores; retirartodaslinhasdeespéciesduplicadasporbloco dámax50=1,55246, comparado comraw1,56426 eordeminvertida1,54716;todos<6,5. Nãoescolhidaslinhasporscore.1000basesROScuradas,min2/mediana164/max194espéciesnãoDOGcomACGT;zerobasesseminformaçãodeoutraespécie. ReciboANO2_duplicate_taxa_sensitivity_20261009_v2/receipt.json. QAinicialfalhouporduplicações,registado,semaprovação;inputseoutputshistóricospreservados. RelatórionovoANO2_REVIEW_AND_DEPENDENCIES_20261009_v3.md. Nãoévalidaçãobiológica nemaprovaçãointegraldaintervenção.

## ANO2 — verificação da anotação UU,09/10/2026
Obtida a anotação NCBI AR106 para GCF_011100685.1, com MD5 oficial conferido e SHA256 guardado. Na janela UU NC_049248.1:6741835–6742835 não há overlap de genes/pseudogenes,exões,CDS ou RNAs. O início5′ anotado mais próximo, incluindo todos os tiposRNA, éNTF3 a101951bp. São fronteiras anotadas, não TSS medidos ou isolamento funcional.
A comparação da lista de risco congelada com símbolos/sinónimos e GeneIDs ROS→UU relacionou2547/2654símbolos;107continuam sem correspondência e não foram transformados em ausência de risco. Entre os genes relacionados, o mais próximo éGAPDH a879877bp;miRNA mais próximo1270499bp. Nenhum hit associado ficou dentro de300kb. Limites e símbolos não resolvidos preservados.
Recibos: INDEPENDENT_LOCI_20261008/ANO2_UU_annotation_crosscheck_20261009/receipt.json, ANO2_UU_risk_miRNA_crosscheck_20261009.json e ANO2_UU_risk_alias_GeneID_crosscheck_20261009.json. Fonte oficial:https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/Canis_lupus_familiaris/106/. Anotação de outroassembly não é confirmação biológica independente. Aprovação integral continua não emitida.


## Fecho da etapa computacional09/10/2026
HAL e UU concluídos; índice UU integral conferido:2197referências/2481983352bases. Duas hipótesesTTTV23 sem hits no modelo definido emROS/UU. Anotação UU revista,sem overlaps;107símbolos de risco permanecem não relacionados,não negativos. ANO2 com suporte aos filtros de locus;aprovação integral não emitida. Stop da seleção de guia:nuclease/variante reais não indicadas. Não é uma exigência de conhecer já o cão experimental. Relatório consolidado:ALL_CANDIDATES_20261005/ANO2_REVIEW_AND_DEPENDENCIES_20261009_v6.md;ledger:INDEPENDENT_LOCI_20261008/ANO2_closure_gate_decision_20261009.json. Nenhum worker ativo desta etapa,contacto externo ou publicação. MemPalace continua emfila local.


## Fecho documentado — 2026-10-12T11:31:25.4215243Z
Relatório ANO2 v7 e decisão v2 concluídos. Cobertura de risco derivada: união 2601/2654, 53 não resolvidos; ROS 2599/55, UU 2592/62. Seis genes locais revistos por GeneID/ortólogos: ANO2/NTF3 sem associação à lista congelada, quatro lncRNAs sem correspondência humana atribuída, mantendo risco desconhecido. Zero novas aprovações integrais. Workers desta etapa terminados. Stop: nuclease/variante real em falta e cobertura de anotação/associação insuficiente para afirmações universais; não é dependência do cão para executar os scans. Automação auditar-todos-os-candidatos-safe-harbor eliminada conforme regra de fecho, sem alegar meta alcançada. MemPalace reavaliado no fecho: SSE404; fila local preservada. Inputs, outputs e contagens históricas intactos; nenhum contacto externo/publicação.
