# Canine safe-harbor candidates — final computational review, 12 October 2026

This is the result of the completed, defined computational stage with the available data. It is not a claim that every possible analysis has been exhausted, that guides are ready for use, or that biological safety has been established. No additional candidate has full audit/intervention approval. The practical locus shortlist is w01 plus ANO2-39708; w11 is not in the advancement panel.

All coordinates below use ROS_Cfam_1.0, BED zero-based, end-exclusive. Windows are 1000 bp. Intervals from other assemblies are not substituted for these coordinates.

| Candidate | ROS interval | Final disposition in this snapshot |
|---|---|---|
| w01 | NC_051811.1:17255706–17256706 | Retained as first priority for testing; no biological validation or complete intervention approval |
| ANO2-39708 | NC_051831.1:39708074–39709074 | Additional provisional pGSH with favourable interrogated locus filters; same broad region as reserve w03, distinct from w01/w11 |
| w11 | NC_051836.1:27041161–27042161 | Not advanced: structural context unresolved, not demonstrated unsafe |
| bg1k_2848 | NC_051840.1:28808728–28809728 | Genomic comparator; blood ATAC P15.21; not a validated T-cell negative control |
| bg1k_0442 | NC_051837.1:10430486–10431486 | Genomic comparator; blood ATAC P7.03; below P10 alone does not establish the combined below-P10 plus 0/76 peak specification |

w01 has 0 annotated repeat bp, exact ROS/UU correspondence, no hit in the queried SV catalogue, and blood ATAC P63.75. Its T-cell accessibility support remains sparse; RRBS observes 20% of the window with mean methylation 60.71%. The nearest annotated non-coding RNA 5′ is 32,785 bp away. Neither low methylation nor functional isolation is demonstrated.

ANO2-39708: full 62,427-bp SHIP-region accessibility P72.70, 20-kb accessibility P65.97; local 1-kb P90.37 is a separate diagnostic. No repeat/regulatory overlap in the queried window sources. Curated maximum 50-bp mean phyloP 1.55246 <6.5. Full ROS/UU sequence searches found zero off-locus hits for two exploratory 23-nt TTTV profiles under the documented substitution/single-bulge/PAM model. These profiles are not approved guides; actual enzyme specificity has not been established. Five local indels and 32 synthetic combinations were evaluated, without asserting phase; 23 offsite catalogue contexts and 40 synthetic combinations did not create model hits. The frozen risk list retains 2654 symbols: 2601 placed in ROS/UU union, 53 unresolved; ROS alone 2599/55, UU alone 2592/62. Four local lncRNAs lack attributed human orthologs in the returned metadata; risk correspondence remains unknown. RRBS covers 15% with mean 48.73%, not whole-window hypomethylation. The full report v7 and its hashed receipts are authoritative for this update.

w11: the catalogued 28.9-Mb simple deletion was not confirmed by coverage, reference support, combined ROS/UU competition, or the reviewed local assemblies/read recruitment. Concordant partial junction sequence is compatible with a local repetitive alteration, but does not reconstruct the complete allele or establish distal linkage. The defined short-read strategy reached a documented resolution limit. The negative conclusion is failure to establish structural suitability, not proof that the deletion is absent or that the locus is unsafe.

Inventory scope: 461 original regions, with 457 retained historical regional exclusions, three regional regulatory conflicts and one insufficient-evidence record. These are regional filter dispositions, not 461 fully assayed local insertion sites. The 64 windows occupy only three broad regions. The 8 October integrated snapshot has 30 unresolved SV-context windows and 34 incomplete reviews, including w01/ANO2. ANO2 follow-ups supersede its earlier pending projection/conservation fields only; the original 64-window snapshot and counts remain unchanged. No other window has been silently promoted.

Expansion: 3000 background windows were screened exploratorily; 133 entered annotation/ATAC triage and 54 a TAD-clear queue. Twenty-three of the 54 have catalogue SV hits; nine combine queried EpiC clearance and no SV hit. Conservation follow-ups completed for bg1k_1059/1266/1900 (0.50206/2.21532/0.46692), but their native gaps do not meet the original convergent 50–75-kb SHIP definition. These are not promoted under the unchanged criteria. A favourable centred extension is not used to replace the native gap.

The earlier regional rankings and counts (including the historical 26 passing regions and the list-overlap argument for independent validation) are superseded by this evidence-bounded review. No genome-wide absence of risk, low methylation, stable CAR expression or T-cell safety is claimed. The remaining enzyme-dependent specificity and annotation-coverage limitations are separate from any subsequent individual or experimental evidence.


## Appendix A — all 64 local windows

These are the immutable 8 October dispositions with an explicit ANO2 follow-up annotation. Functional evidence gaps are not relabelled as demonstrated computational vetoes.

| ROS window ID | Genes | SV hits | 8 October disposition | Pending projection observations | Subsequent interpretation |
|---|---|---:|---|---:|---|
| NC_051811.1_17234956_17235956 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 5 | Original disposition retained; no full audit approval |
| NC_051811.1_17255706_17256706 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 0 | w01 retained for testing, not biologically validated |
| NC_051811.1_17256206_17257206 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 0 | Original disposition retained; no full audit approval |
| NC_051831.1_39713074_39714074 | ANO2/NTF3 | 0 | computational_review_incomplete | 2 | Original disposition retained; no full audit approval |
| NC_051831.1_39729824_39730824 | ANO2/NTF3 | 0 | computational_review_incomplete | 2 | Original disposition retained; no full audit approval |
| NC_051831.1_39730074_39731074 | ANO2/NTF3 | 0 | computational_review_incomplete | 2 | Original disposition retained; no full audit approval |
| NC_051831.1_39730324_39731324 | ANO2/NTF3 | 0 | computational_review_incomplete | 3 | Original disposition retained; no full audit approval |
| NC_051836.1_27039161_27040161 | NPNT/TBCK | 2 | not_advanced_unresolved_SV_context | 24 | Original disposition retained; no full audit approval |
| NC_051836.1_27039411_27040411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 3 | Original disposition retained; no full audit approval |
| NC_051836.1_27040661_27041661 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 4 | Original disposition retained; no full audit approval |
| NC_051836.1_27040911_27041911 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27041161_27042161 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | w11 structural resolution limit documented; not advanced |
| NC_051836.1_27041411_27042411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27041661_27042661 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27050411_27051411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27050661_27051661 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051811.1_17222706_17223706 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 3 | Original disposition retained; no full audit approval |
| NC_051811.1_17222956_17223956 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 2 | Original disposition retained; no full audit approval |
| NC_051811.1_17230206_17231206 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 2 | Original disposition retained; no full audit approval |
| NC_051811.1_17235206_17236206 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 5 | Original disposition retained; no full audit approval |
| NC_051811.1_17235706_17236706 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 1 | Original disposition retained; no full audit approval |
| NC_051811.1_17240706_17241706 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 7 | Original disposition retained; no full audit approval |
| NC_051811.1_17242956_17243956 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 3 | Original disposition retained; no full audit approval |
| NC_051811.1_17245456_17246456 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 0 | Original disposition retained; no full audit approval |
| NC_051811.1_17252456_17253456 | LOC119876429/LOC119872513 | 0 | computational_review_incomplete | 7 | Original disposition retained; no full audit approval |
| NC_051831.1_39682574_39683574 | ANO2/NTF3 | 0 | computational_review_incomplete | 0 | Original disposition retained; no full audit approval |
| NC_051831.1_39685324_39686324 | ANO2/NTF3 | 0 | computational_review_incomplete | 4 | Original disposition retained; no full audit approval |
| NC_051831.1_39685574_39686574 | ANO2/NTF3 | 0 | computational_review_incomplete | 4 | Original disposition retained; no full audit approval |
| NC_051831.1_39686074_39687074 | ANO2/NTF3 | 0 | computational_review_incomplete | 3 | Original disposition retained; no full audit approval |
| NC_051831.1_39686324_39687324 | ANO2/NTF3 | 0 | computational_review_incomplete | 0 | Original disposition retained; no full audit approval |
| NC_051831.1_39686574_39687574 | ANO2/NTF3 | 0 | computational_review_incomplete | 0 | Original disposition retained; no full audit approval |
| NC_051831.1_39707824_39708824 | ANO2/NTF3 | 0 | computational_review_incomplete | 4 | Original disposition retained; no full audit approval |
| NC_051831.1_39708074_39709074 | ANO2/NTF3 | 0 | computational_review_incomplete | 4 | Follow-ups completed in report v7; provisional shortlist; no full approval |
| NC_051831.1_39708324_39709324 | ANO2/NTF3 | 0 | computational_review_incomplete | 2 | Original disposition retained; no full audit approval |
| NC_051831.1_39709574_39710574 | ANO2/NTF3 | 0 | computational_review_incomplete | 0 | Original disposition retained; no full audit approval |
| NC_051831.1_39709824_39710824 | ANO2/NTF3 | 0 | computational_review_incomplete | 0 | Original disposition retained; no full audit approval |
| NC_051831.1_39710074_39711074 | ANO2/NTF3 | 0 | computational_review_incomplete | 14 | Original disposition retained; no full audit approval |
| NC_051831.1_39710324_39711324 | ANO2/NTF3 | 0 | computational_review_incomplete | 14 | Original disposition retained; no full audit approval |
| NC_051831.1_39710574_39711574 | ANO2/NTF3 | 0 | computational_review_incomplete | 14 | Original disposition retained; no full audit approval |
| NC_051831.1_39714074_39715074 | ANO2/NTF3 | 0 | computational_review_incomplete | 1 | Original disposition retained; no full audit approval |
| NC_051831.1_39714824_39715824 | ANO2/NTF3 | 0 | computational_review_incomplete | 10 | Original disposition retained; no full audit approval |
| NC_051831.1_39715074_39716074 | ANO2/NTF3 | 0 | computational_review_incomplete | 8 | Original disposition retained; no full audit approval |
| NC_051831.1_39735074_39736074 | ANO2/NTF3 | 0 | computational_review_incomplete | 16 | Original disposition retained; no full audit approval |
| NC_051836.1_27016161_27017161 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 24 | Original disposition retained; no full audit approval |
| NC_051836.1_27016411_27017411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 24 | Original disposition retained; no full audit approval |
| NC_051836.1_27016661_27017661 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 24 | Original disposition retained; no full audit approval |
| NC_051836.1_27016911_27017911 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 12 | Original disposition retained; no full audit approval |
| NC_051836.1_27017411_27018411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 28 | Original disposition retained; no full audit approval |
| NC_051836.1_27017911_27018911 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 43 | Original disposition retained; no full audit approval |
| NC_051836.1_27018911_27019911 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 30 | Original disposition retained; no full audit approval |
| NC_051836.1_27019911_27020911 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 3 | Original disposition retained; no full audit approval |
| NC_051836.1_27027161_27028161 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 13 | Original disposition retained; no full audit approval |
| NC_051836.1_27027411_27028411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 4 | Original disposition retained; no full audit approval |
| NC_051836.1_27027661_27028661 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 3 | Original disposition retained; no full audit approval |
| NC_051836.1_27033411_27034411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 5 | Original disposition retained; no full audit approval |
| NC_051836.1_27047661_27048661 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27048161_27049161 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27048411_27049411 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27062161_27063161 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 19 | Original disposition retained; no full audit approval |
| NC_051836.1_27062661_27063661 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 3 | Original disposition retained; no full audit approval |
| NC_051836.1_27064161_27065161 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 7 | Original disposition retained; no full audit approval |
| NC_051836.1_27069911_27070911 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 7 | Original disposition retained; no full audit approval |
| NC_051836.1_27070161_27071161 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |
| NC_051836.1_27070204_27071204 | NPNT/TBCK | 1 | not_advanced_unresolved_SV_context | 0 | Original disposition retained; no full audit approval |

## Appendix B — all 461 original regional records

Multiple failed gates may occur in one region; gate counts must not be summed as numbers of distinct regions. Historical records preserved; local-window follow-ups do not globally clear their enclosing region.

| ROS region | Genes | Regional disposition | Failed gates | Historical missing evidence |
|---|---|---|---|---|
| NC_051805.1:5377330–5428711 | LOC111090571/LOC111090572 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:7072137–7132579 | LOC111090579/LOC100685067 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051805.1:12686627–12749345 | LOC111093764/LOC607190 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051805.1:15464961–15523285 | LOC111090610/LOC100688146 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:16284172–16340519 | LOC119869864/PMAIP1 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:23540619–23593624 | LOC119870725/LOC119870724 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:26121059–26180255 | LOC111095273/LOC476208 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:29601219–29667821 | LOC100856045/IL20RA | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:30822469–30875140 | ARFGEF3/LOC111095460 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:32197089–32249732 | LOC119870654/LOC111090654 | excluded_by_legacy_checks | veto_ultraconserved_element; veto_external_regulatory_element; veto_low_atac_accessibility |  |
| NC_051805.1:48451052–48513492 | LOC111095873/RSPH3 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:54446758–54504939 | LOC119869870/PRR18 | excluded_by_legacy_checks | veto_tad_boundary; veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:55288363–55358132 | LOC102153253/LOC111090695 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:55810634–55864353 | LOC119870159/LOC111090696 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:60400908–60462010 | LOC119869937/LOC111089986 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051805.1:65127836–65201393 | TRMT11/LOC100682630 | excluded_by_legacy_checks | veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:68528776–68579924 | LOC111090714/TMEM244 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:71035668–71094862 | LOC111096339/ERCC6L2 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:75730614–75782632 | LOC100684008/RMI1 | excluded_by_legacy_checks | veto_tad_boundary; veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:76492055–76558248 | FRMD3/LOC111096570 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:78699957–78766843 | LOC100685388/LOC476315 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:84616059–84680359 | LOC119869841/LOC100686852 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:85740590–85790734 | LOC119869527/LOC102153097 | excluded_by_legacy_checks | veto_external_regulatory_element |  |
| NC_051805.1:86132349–86196873 | LOC102153343/C1H9orf85 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051805.1:87683331–87744102 | KLF9/SMC5 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:91186888–91251253 | LOC119869582/LOC111096818 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:97204375–97272513 | LOC102156616/LOC111090790 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051805.1:102569417–102626167 | LOC119870888/LOC119871020 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051805.1:118762948–118834530 | LOC100685960/LOC111098307 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051805.1:120003827–120058858 | LOC102152107/DPY19L3 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051806.1:2276949–2331192 | CUL2/LOC111094492 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051806.1:5490645–5562196 | LOC119870094/GAD2 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051806.1:9925163–9998285 | LOC119869979/LOC106560013 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051806.1:13354030–13405141 | LOC111091364/LOC106560015 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051806.1:17055808–17110317 | PTCHD3/ANKRD26 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility |  |
| NC_051806.1:23362298–23415104 | BEND7/LOC119869987 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051806.1:25264511–25315271 | USP6NL/LOC119869990 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051806.1:29613097–29682625 | LOC111091644/LOC111091655 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051806.1:31421264–31473494 | LOC111091701/LOC119870074 | excluded_by_legacy_checks | veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051806.1:37763146–37815322 | LOC102153902/LOC111094200 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051806.1:38503041–38573083 | ARHGAP26/NR3C1 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051806.1:39071176–39129828 | LOC119879638/LOC111092064 | excluded_by_legacy_checks | veto_tad_boundary | repeat_content;conservation;TAD_assignment |
| NC_051806.1:39159547–39230847 | LOC111094604/LOC119870003 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051806.1:55898984–55971473 | PTCD2/ZNF366 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051806.1:61050066–61110203 | LOC119870018/LOC102154386 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051806.1:62723085–62783656 | RPGRIP1L/LOC119870323 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051806.1:63630965–63704659 | TOX3/LOC100687235 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051806.1:77425724–77486897 | LOC102153519/LOC119870436 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051806.1:80612441–80666392 | LOC119870059/KLHDC7A | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051806.1:84152591–84202660 | LOC119870024/PDPN | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051806.1:84303869–84363723 | LRRC38/C2H1orf158 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051807.1:9691905–9747134 | ST8SIA4/LOC106558558 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051807.1:10779807–10837318 | LOC100685686/LOC119871165 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051807.1:21623166–21678171 | LOC119871170/LOC111095350 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:22304729–22367018 | LOC111095151/LOC119871362 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051807.1:25390804–25445201 | LOC119871164/LOC111095157 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:31153307–31211185 | GCNT4/MIR8901 | excluded_by_legacy_checks | veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:32087727–32156058 | TUBGCP5/LOC100688323 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:34846149–34921008 | LOC119871185/LOC106558580 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:35803900–35876584 | UBE3A/LOC100685501 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051807.1:44207725–44265967 | LOC119870785/LOC111095367 | excluded_by_legacy_checks | veto_ultraconserved_element; veto_low_atac_accessibility |  |
| NC_051807.1:50781019–50831091 | AGBL1/LOC102152355 | excluded_by_legacy_checks | veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:56269278–56323122 | EFL1/LOC106558624 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051807.1:60928139–60979861 | LOC119871320/HMX1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051807.1:67665045–67732484 | LOC111095269/LOC111095270 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:72901974–72956918 | RBM47/CHRNA9 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051807.1:73210443–73261055 | LOC119871251/TRNAK-CUU-6 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:74470819–74522415 | LOC106558641/LOC111095389 | excluded_by_legacy_checks | veto_atac_peak; veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051807.1:76113264–76180116 | LOC111095297/LOC102153709 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:76613649–76664485 | ARAP2/LOC100688576 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:77813684–77865855 | LOC111095391/LOC111095392 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051807.1:81859047–81912809 | LOC119872911/LOC119871343 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051807.1:85288823–85342443 | LOC111095403/LOC100686918 | excluded_by_legacy_checks | veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051808.1:1943254–2003893 | LOC102151881/LOC111095689 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:4984288–5034699 | LOC111095431/LOC111095433 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation |
| NC_051808.1:5485501–5546210 | LOC102155069/SLC35F3 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051808.1:6820779–6870891 | LOC119871610/LOC111095435 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:12066504–12128157 | LOC111095700/SLC16A9 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051808.1:13169783–13237184 | CDK1/RHOBTB1 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051808.1:14781335–14832377 | LOC111095463/LOC119871460 | excluded_by_legacy_checks | veto_atac_peak; veto_external_regulatory_element | repeat_content;conservation |
| NC_051808.1:17362530–17419663 | LOC100856460/LOC100686735 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:28998408–29056566 | LOC102155105/LOC111095511 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051808.1:36047237–36120971 | LOC100684047/LOC111095523 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:37518160–37569466 | THOC3/CPLX2 | excluded_by_legacy_checks | veto_risk_gene; veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051808.1:37795829–37851034 | LOC119871736/SFXN1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051808.1:42728958–42791035 | LOC119871511/LOC100687995 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:55513056–55565959 | LOC106558729/LOC119871746 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051808.1:56109927–56183571 | LOC111095578/GRIA1 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051808.1:61210674–61263506 | CSPG4B/LOC119871643 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051808.1:64067111–64119158 | LOC102157144/LOC100687852 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:65949583–66022083 | LOC100687240/MRPS30 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051808.1:76456036–76509373 | DROSHA/CDH6 | excluded_by_legacy_checks | veto_tad_boundary; veto_mirna_nearby | repeat_content;conservation;TAD_assignment |
| NC_051808.1:77033647–77089414 | TMEM125/LOC102156164 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:81940968–82004288 | LOC100685641/LOC479379 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility |  |
| NC_051808.1:83742031–83797839 | CDH12/LOC111095663 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051808.1:87580819–87638502 | ZNF622/LOC106558684 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051809.1:2620428–2685222 | OPCML/NTM | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene | repeat_content;conservation |
| NC_051809.1:6345291–6400601 | LOC111096091/LOC119871991 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051809.1:14259045–14318214 | LOC119871999/LOC111095807 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation |
| NC_051809.1:18909987–18974471 | C5H11orf71/LOC119881489 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051809.1:23237689–23310760 | LOC111095852/C5H11orf87 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051809.1:27846016–27902436 | LOC100687870/LOC111095862 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051809.1:32059172–32124240 | LOC119871841/LOC111095870 | excluded_by_legacy_checks | veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051809.1:47160487–47220638 | LOC119871875/LOC100684598 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051809.1:67812010–67871860 | ZDHHC7/LOC119872032 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051809.1:89203545–89257435 | LOC111096139/LOC119874957 | excluded_by_legacy_checks | veto_tad_boundary | repeat_content;conservation;TAD_assignment |
| NC_051810.1:3695908–3765879 | LOC119881536/LOC119872440 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051810.1:15694571–15754185 | LOC102156691/LOC111096491 | excluded_by_legacy_checks | veto_tad_boundary; veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051810.1:43751242–43809227 | SLC25A24/LOC102151888 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051810.1:46489068–46547626 | LOC111096509/LOC100683041 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051810.1:50658944–50727578 | LOC119872255/LOC119872463 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051810.1:54055016–54116812 | LOC111096377/LOC111096379 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051810.1:54185278–54253504 | LOC119872358/LOC119876387 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051810.1:65303098–65376045 | LOC102151889/LOC111096431 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051810.1:68878658–68947141 | LOC102152406/IFI44 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051811.1:2494758–2545521 | LOC111096593/LOC102152591 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051811.1:11479120–11529669 | LOC111096626/LOC610143 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene | repeat_content;conservation |
| NC_051811.1:15394995–15461327 | LOC106559021/GLUL | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051811.1:17217706–17274511 | LOC119876429/LOC119872513 | regulatory_conflict_requires_local_resolution |  |  |
| NC_051811.1:22837876–22907514 | LOC100688557/LOC111096880 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051811.1:25873242–25943461 | LOC111096660/LOC111096885 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051811.1:35118386–35179633 | LOC111096699/LOC111096700 | excluded_by_legacy_checks | veto_tad_boundary; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051811.1:48020921–48077046 | RIT2/LOC119872716 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051811.1:54309697–54360521 | LOC106559055/GALNT1 | excluded_by_legacy_checks | veto_tad_boundary; veto_mirna_nearby; veto_gene_dense_neighborhood; veto_tad_risk_gene | repeat_content;conservation |
| NC_051811.1:57734073–57792841 | TRAPPC8/LOC111096784 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051811.1:59853346–59916313 | NA/LOC612611 | excluded_by_legacy_checks | veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051811.1:62263653–62326266 | KCTD1/TAF4B | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation |
| NC_051811.1:69954234–70019343 | MYL12B/LOC106559068 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051811.1:72146896–72199972 | LOC111096813/LOC111096814 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051811.1:77245369–77305935 | NA/LOC111096922 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051811.1:80890788–80943582 | MAPK4/LOC100688129 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051812.1:52431–118675 | LOC119872972/LOC119872973 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051812.1:788172–845815 | LOC490585/OR4L1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051812.1:5172690–5239184 | LOC119873153/LOC608162 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051812.1:5534832–5600865 | LOC111097200/LOC119872981 | excluded_by_legacy_checks | veto_ultraconserved_element; veto_low_atac_accessibility |  |
| NC_051812.1:6167980–6227889 | LOC119872983/NOVA1 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051812.1:8640791–8694857 | LOC111096975/LOC102154002 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051812.1:10610167–10663583 | LOC102155565/HEATR5A | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051812.1:12903397–12977974 | NPAS3/EGLN3 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051812.1:21931015–21992476 | LOC100684096/LOC119872816 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility |  |
| NC_051812.1:27735066–27803631 | TMX1/LOC111097216 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051812.1:28781785–28839966 | LOC111097218/TXNDC16 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051812.1:36266983–36333691 | SLC38A6/LOC111097065 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation |
| NC_051812.1:54071776–54128751 | LOC102153790/LOC119873009 | excluded_by_legacy_checks | veto_tad_boundary; veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051812.1:55781421–55851186 | LOC111097231/LOC106559104 | excluded_by_legacy_checks | veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051812.1:58732540–58805036 | LOC119873149/LOC119873016 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility |  |
| NC_051812.1:61680835–61737486 | CALM1/LOC111097143 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051812.1:67788194–67843405 | LOC119873039/LOC111097248 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene | repeat_content;conservation |
| NC_051812.1:70725627–70776316 | LOC106559203/LOC119881569 | excluded_by_legacy_checks | veto_atac_peak; veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051812.1:71254791–71309046 | TNFAIP2/LOC480441 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051813.1:8313493–8376479 | LOC102156472/LOC102156400 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051813.1:8406782–8479948 | LOC119873409/SOX9 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051813.1:15175682–15238845 | LOC102154199/CACNG5 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051813.1:33214836–33269418 | MSI2/CCDC182 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene | repeat_content;conservation |
| NC_051813.1:35794767–35867117 | BRIP1/LOC100855486 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051813.1:42596793–42666787 | LOC102151497/LOC111097438 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051813.1:44260433–44316601 | LOC119873524/LOC119881587 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051813.1:53233104–53304110 | MED27/LOC111097551 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051813.1:55289222–55342350 | LOC102151894/LOC111097492 | excluded_by_legacy_checks | veto_atac_peak; veto_mirna_nearby; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051813.1:60390416–60460489 | LOC119873444/LOC119873377 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:6647061–6715959 | LOC119881620/LOC119881618 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051814.1:16586697–16654416 | LOC100685032/LOC119873751 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051814.1:22082493–22156729 | LOC119873757/LOC102153700 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:24918877–24986657 | LOC119873626/LOC119873758 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:25211889–25276787 | SLC25A17/MCHR1 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:29447517–29507386 | LOC102152977/MCM5 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:33387381–33437530 | LOC106559418/LOC106559383 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_tad_risk_gene | repeat_content;conservation |
| NC_051814.1:34678177–34747422 | LOC102152060/TXNRD1 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:35158873–35212398 | SEPTIN10/LOC106559442 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051814.1:46275184–46345564 | TRNAI-UAU-2/LOC119873662 | excluded_by_legacy_checks | veto_tad_boundary; veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051814.1:46556671–46612324 | LOC106559407/LOC119873664 | excluded_by_legacy_checks | veto_atac_peak; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:55173498–55233905 | NA/LOC102151272 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051814.1:57054512–57104779 | LOC119873779/CLHC1 | excluded_by_legacy_checks | veto_atac_peak; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051814.1:57815583–57880350 | LOC119873681/LOC111097852 | excluded_by_legacy_checks | veto_mirna_nearby; veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051814.1:61932847–61989787 | LOC102155693/LOC119873691 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene | repeat_content;conservation |
| NC_051814.1:63844443–63896703 | OTX1/LOC100686739 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051814.1:66943912–67008116 | MEIS1/LOC111097797 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051814.1:68117840–68179368 | LOC100683563/C1D | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051815.1:5132068–5205139 | LOC102157398/LOC119873853 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051815.1:7585907–7653876 | LOC100687242/LOC106559451 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility |  |
| NC_051815.1:15804115–15863613 | LOC111097925/LOC111098006 | excluded_by_legacy_checks | veto_ultraconserved_element; veto_external_regulatory_element; veto_low_atac_accessibility |  |
| NC_051815.1:33774377–33832861 | LOC100856273/LOC100856119 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051815.1:54176384–54240721 | MELK/LOC119874128 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051815.1:56181459–56247618 | LOC111098102/LOC119874060 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051815.1:67619694–67676552 | LOC119873978/LOC111097883 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051815.1:72532281–72586362 | LOC106559529/LOC119873989 | excluded_by_legacy_checks | veto_tad_boundary; veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051815.1:74432972–74494463 | LOC111097897/LOC119874053 | excluded_by_legacy_checks | veto_ultraconserved_element; veto_external_regulatory_element |  |
| NC_051816.1:729204–780179 | LOC100684529/LOC119874351 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051816.1:888522–959052 | MUC21/LOC119876866 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051816.1:2021125–2084815 | LOC119881685/TSBP1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:5010995–5063926 | LOC608284/FKBP5 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:10056117–10121857 | LOC111098152/LOC119874267 | excluded_by_legacy_checks | veto_risk_gene_radius | repeat_content;conservation |
| NC_051816.1:13571280–13642018 | LOC111098270/SUPT3H | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius | repeat_content;conservation |
| NC_051816.1:16360664–16416333 | OPN5/LOC119876862 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:25815845–25884245 | LOC111098165/LOC111098166 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:40263604–40321877 | LOC100855752/LOC111098325 | excluded_by_legacy_checks | veto_high_repeat_content; veto_low_atac_accessibility |  |
| NC_051816.1:41991512–42044662 | BCKDHB/LOC102155357 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:47327902–47386565 | HTR1E/CGA | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:48272208–48342011 | SPACA1/CNR1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:48393708–48452461 | LOC111098352/LOC119874280 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:53806181–53865232 | LOC100684677/LOC100687875 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:64340623–64411071 | PRDM1/ATG5 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051816.1:67608124–67672326 | GPR6/LOC111098188 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051816.1:70866191–70918704 | MARCKS/LOC106559557 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051816.1:72343871–72396378 | LOC119874236/LOC111098409 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051816.1:73372682–73425689 | LOC111098416/PSMB1 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051817.1:670319–729017 | LOC111098485/LOC111098432 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051817.1:2762900–2833447 | LOC119876950/LOC119864717 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051817.1:8064988–8132108 | LOC111098516/LOC119865560 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051817.1:8580007–8642797 | LOC111098518/LOC111098519 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051817.1:10033333–10096852 | MIR8836/NUDCD1 | excluded_by_legacy_checks | veto_mirna_nearby; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051817.1:12626538–12686868 | LOC106559613/CSMD3 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility |  |
| NC_051817.1:16507341–16565029 | UTP23/LOC100687037 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051817.1:24350676–24411685 | LOC119865453/LOC102151612 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051817.1:28001742–28068918 | LOC111098451/ADCY8 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051817.1:29429036–29482883 | LOC102153335/LRRC6 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051817.1:42409351–42467972 | LOC119865519/LOC119865447 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051817.1:50196685–50266401 | LOC100683751/LOC100683821 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051818.1:5712360–5767938 | LOC102155278/LOC119867007 | excluded_by_legacy_checks | veto_tad_boundary; veto_mirna_nearby; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051818.1:9869905–9920386 | LOC119865745/LOC100688839 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051818.1:14495427–14554191 | LOC119874506/TEX47 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051818.1:19245636–19304857 | GNG11/BET1 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051818.1:43701564–43775013 | ADCYAP1R1/LOC102154587 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051818.1:47079968–47131672 | LOC119867666/LOC111098822 | excluded_by_legacy_checks | veto_risk_gene_radius | repeat_content;conservation |
| NC_051818.1:60626704–60697896 | LOC119866869/SLC13A1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051819.1:3753627–3825487 | LOC119869470/LOC119881733 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051819.1:4205082–4273748 | LOC111098951/LOC102151429 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051819.1:20314140–20377667 | LOC119869464/LOC119869379 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051819.1:23412529–23465999 | MYF5/LIN7A | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051819.1:27740518–27794146 | LOC102154855/LOC100683898 | excluded_by_legacy_checks | veto_mirna_nearby; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051819.1:27864256–27917940 | LOC111090020/RASSF9 | excluded_by_legacy_checks | veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051819.1:29615719–29672503 | LOC111098932/LOC111098894 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051819.1:30913325–30967876 | LOC119869342/DUSP6 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051819.1:32255898–32317043 | LOC111090028/CCER1 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051819.1:35221213–35288244 | LOC119869290/TMCC3 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051819.1:50452875–50525877 | LOC119869387/LOC111098915 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051819.1:56330968–56393135 | LOC119867837/LOC119869360 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051820.1:1623305–1690052 | PDIA4/LOC119869648 | excluded_by_legacy_checks | veto_atac_peak; veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051820.1:12347409–12407429 | LUZP6/LOC102156014 | excluded_by_legacy_checks | veto_atac_peak; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051820.1:22590144–22660334 | DYNC2I1/VIPR2 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051820.1:22976673–23048542 | LOC111090199/LOC100682550 | insufficient_evidence |  | conservation |
| NC_051820.1:35121001–35190439 | LOC100686540/WRN | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051820.1:37178242–37231175 | LOC119869555/LOC102155995 | excluded_by_legacy_checks | veto_tad_boundary; veto_mirna_nearby; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051820.1:37416979–37473636 | PPP1R3B/ERI1 | excluded_by_legacy_checks | veto_mirna_nearby; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051820.1:38028250–38082724 | LOC111090314/LOC111090316 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051820.1:46125770–46187174 | LOC119869631/LOC100682827 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051820.1:46199665–46254098 | LOC119869583/LOC119869626 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051820.1:49468828–49528532 | LOC119869512/LOC102156721 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051820.1:51445768–51508808 | LOC119874730/LOC119874609 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051820.1:54629014–54692819 | AGA/NEIL3 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051821.1:4818014–4875944 | LOC111090412/LOC102151619 | excluded_by_legacy_checks | veto_external_regulatory_element |  |
| NC_051821.1:11891971–11964254 | LOC102157135/LOC111090504 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051821.1:12278057–12338750 | LOC102157261/CYRIA | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051821.1:15968485–16030074 | LOC111090394/APOB | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051821.1:17141980–17213143 | LOC111090524/LOC119881779 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051821.1:19088682–19138909 | LOC111090423/LOC111090532 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051821.1:25527391–25580014 | LOC119869680/SRD5A2 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051821.1:34751027–34802687 | LOC102151793/KCNG3 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051821.1:51685511–51738942 | LOC119869763/LOC111090442 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051821.1:52870289–52932765 | LOC111090443/LOC119864061 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius | repeat_content;conservation |
| NC_051821.1:52940905–52993050 | LOC102154674/TRIM33 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius | repeat_content;conservation;TAD_assignment |
| NC_051821.1:53716583–53781749 | LOC100688064/NGF | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051821.1:55900409–55960106 | LOC111090684/LOC100688781 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051821.1:57904516–57954590 | LOC102154106/LOC111090694 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051821.1:62647649–62699011 | LOC119864001/LOC119878526 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051821.1:62708801–62760445 | LOC119878525/LOC119877138 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051822.1:7159763–7220802 | LOC111090865/LOC102151924 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051822.1:11006057–11065538 | VPS41/LOC119864320 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051822.1:12195480–12258671 | LOC100683828/EPDR1 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051822.1:18243340–18307563 | PHTF2/MAGI2 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051822.1:23806419–23872683 | LOC111090917/SEMA3A | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051822.1:30797965–30858164 | LOC111091010/LOC100682696 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051822.1:36175607–36236940 | LOC111090814/LOC119864218 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene | repeat_content;conservation |
| NC_051822.1:37530345–37599747 | MPPED2/ARL14EP | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051822.1:37625161–37675210 | LOC111090955/OR10V5 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051822.1:55169332–55241058 | AHNAK/SCGB1A1 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051822.1:55656583–55708852 | LOC100856625/LOC111091015 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051823.1:15192889–15259450 | LOC119864516/FAT4 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051823.1:22718008–22788982 | NA/LOC111091051 | excluded_by_legacy_checks | veto_tad_boundary; veto_atac_peak | repeat_content;conservation |
| NC_051823.1:32708863–32761746 | NA/LOC100683928 | excluded_by_legacy_checks | veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051823.1:40699370–40753962 | LOC119864393/LOC119864525 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051823.1:42150441–42214374 | LOC119864487/LOC119864518 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051823.1:43660794–43734482 | LOC119864395/LOC119864396 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051823.1:46607767–46663039 | KYNU/LOC111091185 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051823.1:50766576–50830568 | LOC111091082/LOC111091197 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051823.1:51746959–51816629 | LOC111091200/LOC119864405 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051823.1:53515801–53574039 | LOC119864407/LOC119864408 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051823.1:53669125–53730654 | LOC100687327/LOC119864410 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051824.1:4223236–4289981 | LOC100683648/CHCHD4 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051824.1:11749774–11807289 | LOC100685632/LOC111091252 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051824.1:12427622–12487484 | LOC111091257/EDEM1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051824.1:13377340–13441692 | LOC476550/LOC100686411 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051824.1:17040502–17103820 | LOC102154165/LOC119864650 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051824.1:19442546–19500025 | SHQ1/LOC102154862 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element | repeat_content;conservation |
| NC_051824.1:21700981–21759122 | MDFIC2/LOC119864576 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051824.1:24572625–24634330 | LOC102157163/LOC111091281 | excluded_by_legacy_checks | veto_tad_boundary; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051824.1:32811725–32865996 | LOC111091305/SLMAP | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051825.1:1063011–1131564 | LOC111091473/LOC606831 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051825.1:11676400–11742674 | GRM5/LOC119864956 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051825.1:15952393–16002786 | DLG2/LOC100686619 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051825.1:38848629–38901890 | LOC102152988/CRSP-2 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051825.1:48906661–48959809 | BBOX1/LOC111091586 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051825.1:51990876–52054339 | LOC119864909/LOC119864867 | excluded_by_legacy_checks | veto_mirna_nearby; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051826.1:1416239–1480570 | LOC111091743/LOC111091793 | excluded_by_legacy_checks | veto_external_regulatory_element |  |
| NC_051826.1:5334971–5405776 | ZC3H13/LOC111091811 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051826.1:10143251–10217071 | LOC111091828/LOC100685259 | excluded_by_legacy_checks | veto_mirna_nearby | repeat_content;conservation |
| NC_051826.1:14184578–14251277 | PCDH17/LOC119865129 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051826.1:16687132–16738584 | TDRD3/LOC100686550 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051826.1:19596764–19663367 | LOC100687588/LOC100687653 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051826.1:27517646–27576124 | KLF5/LOC100682890 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051826.1:29235878–29299159 | LOC100684181/LOC111091654 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051826.1:33638670–33695479 | LOC106557558/SPRY2 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051826.1:36431283–36491812 | LOC100684683/LOC111091757 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051826.1:45865441–45937672 | GPR180/SOX21 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051826.1:55453706–55512794 | LOC100683697/LOC111091706 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051826.1:55859005–55918918 | LOC119865267/LOC119865154 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051826.1:59581650–59643922 | ING1/LOC119865109 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051827.1:1109045–1172665 | SYNDIG1/LOC100688361 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051827.1:2619236–2673562 | LOC100687334/LOC111091893 | excluded_by_legacy_checks | veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051827.1:3227875–3294390 | LOC111091943/LOC119865296 | excluded_by_legacy_checks | veto_atac_peak; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051827.1:6808405–6874087 | STAC/LOC100685633 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051827.1:10487084–10553541 | LOC111091857/LOC106557623 | excluded_by_legacy_checks | veto_tad_risk_gene | repeat_content;conservation |
| NC_051827.1:23551824–23607775 | LOC111091912/LOC119865327 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051827.1:29008575–29073911 | NUDT16/MRPL3 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051827.1:34198197–34269938 | LOC111091952/LOC119865310 | excluded_by_legacy_checks | veto_ultraconserved_element; veto_low_atac_accessibility |  |
| NC_051827.1:34370653–34427811 | LOC111091954/LOC119865330 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051827.1:39231375–39301877 | CHST2/LOC100688497 | excluded_by_legacy_checks | veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051827.1:45474517–45526737 | LOC100855515/LOC119865420 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051828.1:11487994–11545704 | LOC119865622/LOC111092144 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051828.1:16606256–16661497 | GPCPD1/LOC119865589 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051828.1:19302591–19370767 | LOC111092201/STK35 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051828.1:26717079–26778334 | SRC/LOC111092200 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051828.1:36042329–36111036 | LOC119865683/PREX1 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051828.1:36665388–36721820 | LOC111092109/KCNB1 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051828.1:37513357–37574621 | LOC119865562/LOC102152659 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation |
| NC_051828.1:40064815–40121938 | TSHZ2/ZNF217 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051828.1:42475068–42526723 | LOC611771/LOC119865771 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051828.1:42745535–42819449 | LOC119865786/LOC119865699 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051829.1:326595–381037 | THSD1/LOC100684545 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051829.1:804447–863659 | FOXO1/LOC102153908 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051829.1:9416764–9476935 | LOC119865997/LOC119865829 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051829.1:12228903–12283288 | LOC106557742/LOC102153058 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051829.1:16583402–16638128 | LOC119865834/LOC106557750 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051829.1:20424770–20475505 | LOC111092453/NEK1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051829.1:34024080–34083975 | NKX3-1/LOC111092402 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene | repeat_content;conservation |
| NC_051829.1:40554664–40610407 | AGFG1/LOC486150 | excluded_by_legacy_checks | veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051829.1:45854728–45924770 | LOC119866032/LOC119865861 | excluded_by_legacy_checks | veto_atac_peak; veto_external_regulatory_element | repeat_content;conservation |
| NC_051829.1:46659216–46733469 | LOC111092437/LOC119865863 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051830.1:3084750–3142913 | SLC15A4/TMEM132C | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood | repeat_content;conservation;TAD_assignment |
| NC_051830.1:4496506–4546999 | LOC119866086/LOC477441 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051830.1:12974593–13028567 | LOC119866092/MED13L | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element | repeat_content;conservation |
| NC_051830.1:17193803–17254161 | LOC119866085/LOC119866121 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051830.1:20964282–21018696 | LOC119866095/LOC119866096 | excluded_by_legacy_checks | veto_tad_risk_gene | repeat_content;conservation |
| NC_051830.1:37136455–37205993 | LOC106557850/LOC102154823 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051831.1:11296029–11362267 | LOC106557871/LOC111092904 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051831.1:12006551–12065044 | GXYLT1/LOC111092835 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051831.1:13274723–13345552 | NA/MUC19 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051831.1:17823426–17886691 | IPO8/LOC100685338 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051831.1:32573520–32640058 | LOC102155441/LOC100685365 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element | repeat_content;conservation |
| NC_051831.1:39677324–39739751 | ANO2/NTF3 | regulatory_conflict_requires_local_resolution |  |  |
| NC_051831.1:40299441–40358999 | LOC119866433/GALNT8 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051832.1:2510899–2563279 | NA/OR6D5 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051832.1:6155210–6215233 | LOC119866635/LOC106557967 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051832.1:24424569–24481957 | LOC100855860/LOC111093048 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051832.1:25443005–25493730 | ADRB1/CCDC186 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051832.1:34214445–34267922 | LOC111093045/LOC119866531 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation |
| NC_051833.1:357464–431756 | UBE2V2/LOC100856084 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051833.1:1669770–1738756 | LOC100684688/LOC119866769 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051833.1:20829395–20898312 | LOC119866675/MSC | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051833.1:21906082–21964104 | LOC111093217/LOC119866789 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051833.1:23112498–23172320 | LOC111093212/LOC111093135 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051833.1:23509733–23568950 | CRISPLD1/LOC100688946 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051833.1:28263136–28333785 | LOC106558008/PAG1 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051833.1:28623707–28681799 | FABP5/LOC119866803 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051833.1:38580004–38634159 | LOC100685448/LOC119866804 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051834.1:2536913–2592033 | LOC111093357/GJD2 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051834.1:4478505–4544211 | CDIN1/MEIS2 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051834.1:11061036–11119847 | LOC119866880/LOC111093393 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051834.1:24994850–25055410 | FOXB1/LOC119866903 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051834.1:32644888–32717232 | LOC111093320/LOC100684041 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood | repeat_content;conservation;TAD_assignment |
| NC_051834.1:34669412–34723012 | LOC102154511/LARP6 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051834.1:38942658–39000282 | LINGO1/LOC119866844 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051835.1:955347–1016678 | LOC100688275/LOC102152193 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051835.1:13106261–13156978 | LOC119867042/LOC111093424 | excluded_by_legacy_checks | veto_mirna_nearby; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051835.1:18249632–18323889 | LOC100685420/LOC100685577 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051835.1:23436373–23508916 | LOC111093569/LOC119867012 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility |  |
| NC_051835.1:33017236–33069006 | LOC102152343/LOC111093484 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051835.1:33785714–33837005 | LOC119867070/LOC111093482 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood | repeat_content;conservation;TAD_assignment |
| NC_051835.1:36406001–36461046 | SLC37A1/LOC119867027 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051836.1:3698661–3752870 | LOC111093630/GK2 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051836.1:11276866–11345764 | MEPE/LOC111093643 | excluded_by_legacy_checks | veto_atac_peak; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051836.1:27013911–27071204 | NPNT/TBCK | regulatory_conflict_requires_local_resolution |  |  |
| NC_051837.1:574972–634483 | LOC111093784/LOC111093791 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051837.1:7790832–7842018 | LOC111093705/IMPG2 | excluded_by_legacy_checks | veto_low_atac_accessibility |  |
| NC_051837.1:8568536–8621477 | LOC119867284/LOC111093752 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051837.1:11353405–11418044 | ALCAM/CBLB | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051837.1:20017414–20079909 | GAP43/LSAMP | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051837.1:31269207–31329671 | LOC106558167/FAM43A | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051838.1:7080793–7147581 | NSUN2/UBE2QL1 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051838.1:10054548–10126362 | LOC111093922/LOC119867475 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051838.1:19004008–19067846 | LOC119867441/CRYGS | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051838.1:19321861–19376696 | ADIPOQ/LOC119867406 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051838.1:33370196–33420674 | MIR551B/LOC119867412 | excluded_by_legacy_checks | veto_mirna_nearby; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051838.1:33641835–33705242 | LOC119867455/LOC102151967 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051838.1:34380069–34430694 | LOC102152454/LOC111093889 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051838.1:36725342–36775647 | FNDC3B/GHSR | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051838.1:40621835–40683716 | LOC111093930/TBL1XR1 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element | repeat_content;conservation;TAD_assignment |
| NC_051839.1:1643004–1708560 | LOC119867596/LOC119867597 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051839.1:1805405–1856208 | FOXF2/LOC100684996 | excluded_by_legacy_checks | veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051839.1:11518999–11569382 | LOC111093992/LOC111094013 | excluded_by_legacy_checks | veto_atac_peak; veto_external_regulatory_element | repeat_content;conservation |
| NC_051839.1:18052245–18123443 | LOC111093969/LOC102156827 | excluded_by_legacy_checks | veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051839.1:18480986–18533175 | LOC111094025/LOC111094007 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051840.1:1166863–1219794 | GALNT13/LOC106558288 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051840.1:7616013–7676372 | LOC111094266/GCG | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051840.1:8466112–8532234 | LOC100688550/LOC119867785 | excluded_by_legacy_checks | veto_low_atac_accessibility | repeat_content;conservation |
| NC_051840.1:14845639–14901149 | UBR3/LOC111094278 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051840.1:16124607–16177018 | CYBRD1/LOC111094279 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood | repeat_content;conservation |
| NC_051840.1:16788107–16853124 | LOC111094247/LOC119867830 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051840.1:24054229–24128740 | LOC119867869/LOC111094258 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051840.1:24296513–24361259 | LOC111094259/LOC119867786 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051840.1:26520111–26586822 | LOC111094249/LOC102155713 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051841.1:2245776–2303016 | LOC119867916/TMEFF2 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element | repeat_content;conservation |
| NC_051841.1:7533022–7602083 | LOC100686361/LOC100686866 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility |  |
| NC_051841.1:12473128–12532864 | CD28/LOC100856338 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene | repeat_content;conservation |
| NC_051841.1:27629161–27690392 | LOC119867980/EPHA4 | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051841.1:28495545–28545686 | LOC106558331/LOC119867969 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene | repeat_content;conservation |
| NC_051841.1:29390284–29459285 | LOC111094442/SCG2 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051841.1:29837431–29901354 | LOC111094360/LOC119867910 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051842.1:3541738–3613722 | LOC111094553/LOC111094552 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051842.1:6148634–6217613 | LOC119868013/RGS2 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051842.1:8585823–8643867 | BRINP3/LOC119868094 | excluded_by_legacy_checks | veto_tad_boundary; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051842.1:16412131–16486695 | LOC111094583/LOC111094562 | excluded_by_legacy_checks | veto_external_regulatory_element | repeat_content;conservation |
| NC_051842.1:19592201–19658874 | LOC111094531/NUF2 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051842.1:22834825–22886249 | OR10J2/ACKR1 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051843.1:6868917–6934242 | CLCN4/LOC119868598 | excluded_by_legacy_checks | veto_tad_boundary; veto_external_regulatory_element | repeat_content;conservation |
| NC_051843.1:8354252–8409529 | LOC106558419/LOC119868687 | excluded_by_legacy_checks | veto_low_mappability; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:10578732–10643009 | LOC480844/FANCB | excluded_by_legacy_checks | veto_external_regulatory_element; veto_low_atac_accessibility |  |
| NC_051843.1:15197385–15269910 | LOC119868379/LOC119868161 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:26196588–26269373 | TRNAE-CUC-22/LOC119868188 | excluded_by_legacy_checks | veto_tad_boundary; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051843.1:32334636–32387595 | LOC612281/LOC111094893 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation;TAD_assignment |
| NC_051843.1:34617972–34669640 | LOC119868190/BCOR | excluded_by_legacy_checks | veto_risk_gene; veto_low_mappability; veto_risk_gene_radius; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:35941631–35992808 | LOC111094865/CASK | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:50115742–50188659 | LOC119868203/WAS | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:51534278–51599514 | LOC100686592/LOC491922 | excluded_by_legacy_checks | veto_mirna_nearby; veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:56733328–56786162 | TAF1/LOC100688354 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:59813343–59874197 | UPRT/ZDHHC15 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_high_repeat_content; veto_low_atac_accessibility |  |
| NC_051843.1:69851587–69907041 | PABPC5/LOC111094733 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_high_repeat_content; veto_low_atac_accessibility |  |
| NC_051843.1:74667098–74718660 | LOC119868578/LOC119868214 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:77845859–77904191 | BHLHB9/LOC119868603 | excluded_by_legacy_checks | veto_mirna_nearby; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:79168554–79226339 | FAM199X/ESX1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility |  |
| NC_051843.1:82817447–82882230 | PRPS1/TSC22D3 | excluded_by_legacy_checks | veto_atac_peak; veto_gene_dense_neighborhood; veto_external_regulatory_element | repeat_content;conservation |
| NC_051843.1:90957857–91025361 | LOC119868706/LOC111094871 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:97127995–97179954 | XIAP/LOC119868511 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_gene_dense_neighborhood; veto_tad_risk_gene; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:99485997–99542726 | LOC492115/LOC609540 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:104147730–104215864 | LOC119868710/IGSF1 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility |  |
| NC_051843.1:105250253–105306545 | LOC119868568/LOC111094836 | excluded_by_legacy_checks | veto_tad_boundary; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:105886803–105957615 | LOC111094882/GPC4 | excluded_by_legacy_checks | veto_risk_gene_radius; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:116212579–116266672 | LOC100686375/LOC100856426 | excluded_by_legacy_checks | veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:122201824–122251894 | LOC100686136/LOC492221 | excluded_by_legacy_checks | veto_mirna_nearby; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:122735338–122788464 | LOC612447/LOC119868662 | excluded_by_legacy_checks | veto_mirna_nearby; veto_gene_dense_neighborhood; veto_low_atac_accessibility | repeat_content;conservation |
| NC_051843.1:122883327–122939331 | LOC119868539/CETN2 | excluded_by_legacy_checks | veto_mirna_nearby; veto_gene_dense_neighborhood; veto_external_regulatory_element; veto_low_atac_accessibility | repeat_content;conservation |

## Appendix C — all 54 independent exploratory queue windows

Original queue flags are screening results, not final approvals. The three completed conservation follow-ups are discussed above; favourable subtests do not override SHIP geometry.

| ID | ROS interval | SV hits | EpiC clear/no SV with roundtrip | Original approved flag |
|---|---|---:|---|---|
| bg1k_0028 | NC_051833.1:11531417–11532417 | 1 | False | False |
| bg1k_0032 | NC_051822.1:31499321–31500321 | 1 | False | False |
| bg1k_0182 | NC_051809.1:28289785–28290785 | 0 | False | False |
| bg1k_0186 | NC_051810.1:34191173–34192173 | 1 | False | False |
| bg1k_0209 | NC_051806.1:71237525–71238525 | 0 | False | False |
| bg1k_0266 | NC_051808.1:64832915–64833915 | 0 | False | False |
| bg1k_0272 | NC_051830.1:35147473–35148473 | 1 | False | False |
| bg1k_0285 | NC_051810.1:54092951–54093951 | 1 | False | False |
| bg1k_0313 | NC_051828.1:41753808–41754808 | 0 | False | False |
| bg1k_0334 | NC_051825.1:16855903–16856903 | 0 | False | False |
| bg1k_0403 | NC_051810.1:34396997–34397997 | 1 | False | False |
| bg1k_0455 | NC_051805.1:121281018–121282018 | 0 | False | False |
| bg1k_0529 | NC_051805.1:9698901–9699901 | 2 | False | False |
| bg1k_0626 | NC_051807.1:10280149–10281149 | 0 | False | False |
| bg1k_0631 | NC_051807.1:34514227–34515227 | 0 | False | False |
| bg1k_0809 | NC_051826.1:4030104–4031104 | 1 | False | False |
| bg1k_0836 | NC_051806.1:32999269–33000269 | 0 | True | False |
| bg1k_0848 | NC_051839.1:13978557–13979557 | 1 | False | False |
| bg1k_0861 | NC_051809.1:26893945–26894945 | 0 | True | False |
| bg1k_0949 | NC_051826.1:56418779–56419779 | 1 | False | False |
| bg1k_0984 | NC_051822.1:31564724–31565724 | 1 | False | False |
| bg1k_1043 | NC_051807.1:6162258–6163258 | 0 | False | False |
| bg1k_1059 | NC_051831.1:19034940–19035940 | 0 | True | False |
| bg1k_1186 | NC_051827.1:4735013–4736013 | 1 | False | False |
| bg1k_1206 | NC_051817.1:24272763–24273763 | 0 | False | False |
| bg1k_1216 | NC_051840.1:23998261–23999261 | 0 | True | False |
| bg1k_1217 | NC_051831.1:18142232–18143232 | 1 | False | False |
| bg1k_1231 | NC_051828.1:41184704–41185704 | 0 | False | False |
| bg1k_1266 | NC_051807.1:2091307–2092307 | 0 | True | False |
| bg1k_1267 | NC_051816.1:9525800–9526800 | 0 | False | False |
| bg1k_1317 | NC_051822.1:30325702–30326702 | 1 | False | False |
| bg1k_1367 | NC_051843.1:96269013–96270013 | 1 | False | False |
| bg1k_1386 | NC_051837.1:3475382–3476382 | 0 | False | False |
| bg1k_1526 | NC_051809.1:11066492–11067492 | 1 | False | False |
| bg1k_1654 | NC_051806.1:45531845–45532845 | 0 | True | False |
| bg1k_1786 | NC_051824.1:10227253–10228253 | 1 | False | False |
| bg1k_1796 | NC_051826.1:8423204–8424204 | 2 | False | False |
| bg1k_1865 | NC_051812.1:32268116–32269116 | 1 | False | False |
| bg1k_1869 | NC_051815.1:49703204–49704204 | 0 | False | False |
| bg1k_1900 | NC_051814.1:43352999–43353999 | 0 | True | False |
| bg1k_1922 | NC_051823.1:34462065–34463065 | 0 | False | False |
| bg1k_2150 | NC_051807.1:84418156–84419156 | 0 | True | False |
| bg1k_2163 | NC_051817.1:28742272–28743272 | 0 | False | False |
| bg1k_2176 | NC_051829.1:32973602–32974602 | 0 | False | False |
| bg1k_2275 | NC_051823.1:28387641–28388641 | 0 | False | False |
| bg1k_2292 | NC_051805.1:47450307–47451307 | 2 | False | False |
| bg1k_2340 | NC_051816.1:43724382–43725382 | 0 | True | False |
| bg1k_2382 | NC_051827.1:26596631–26597631 | 0 | False | False |
| bg1k_2403 | NC_051821.1:27458567–27459567 | 0 | False | False |
| bg1k_2427 | NC_051823.1:38132810–38133810 | 0 | False | False |
| bg1k_2448 | NC_051812.1:44637049–44638049 | 2 | False | False |
| bg1k_2544 | NC_051805.1:1486678–1487678 | 0 | False | False |
| bg1k_2676 | NC_051821.1:48836825–48837825 | 1 | False | False |
| bg1k_2838 | NC_051830.1:13454672–13455672 | 2 | False | False |
