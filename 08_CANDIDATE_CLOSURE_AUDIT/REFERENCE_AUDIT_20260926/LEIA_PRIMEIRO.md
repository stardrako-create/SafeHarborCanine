# Auditoria de referência — 26/09/2026

O cabeçalho do VCF não contém declaração reference/assembly. Os comprimentos dos 38 autossomas coincidem com UU e nenhum coincide com ROS. Chr32: catálogo/UU 40.225.481 bp; ROS 39.023.732 bp. Isto verifica nomenclatura/comprimentos, não identidade integral da sequência.

A fonte primária Dog10K identifica UU_Cfam_GSD_1.0 (GCF_011100685.1), suplementada com três contigs Y de ROS. Fonte: https://pmc.ncbi.nlm.nih.gov/articles/PMC10426128/ ; nota de distribuição: https://kiddlabshare.med.umich.edu/dog10K/Dog10K-data-release-note.pdf . O FASTA local UU não foi aqui demonstrado idêntico à referência suplementada completa; a concordância autossómica não resolve diferenças em Y/contigs adicionais.

w11 UU NC_049253.1:13078664–13079664 alinha 1000M, NM0, em ROS NC_051836.1:27041161–27042161, orientação inversa (BED0).

Os flancos de10kb dos dois limites alinham predominantemente ao chr32 ROS em orientação inversa, um próximo de35,34Mb e outro de6,41Mb. Há clipping e alinhamento suplementar no flanco direito; não se fez liftOver exato de todos os breakpoints.

A junção sintética de2kb não tem alinhamento contínuo integral reportado em ROS: os dois segmentos de cerca de1kb alinham separadamente nas regiões distantes acima. Portanto, esta comparação não apoia a explicação simples de que a junção seja apenas a organização contínua normal de ROS. BWA é heurístico: ausência de alinhamento reportado não prova inexistência de todas as alternativas.

Usar UU para investigar a chamada do catálogo foi coerente com a sua referência. ROS permanece a referência principal do candidato. A ligação explícita entre as análises estava incompleta e foi agora documentada. A variante estrutural continua sem resolução; este resultado não estabelece deleção nem genótipo.

Ficheiros: assembly_lengths.json, queries.fa, queries_ROS.sam, bwa.log e review.json. Inputs anteriores preservados. Não houve nova execução WGS nem reativação da automação.
