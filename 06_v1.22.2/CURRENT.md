# Safe Harbor CAR-T — v1.22.2: correspondência LOC e detalhe de mapeabilidade NPNT

15 de setembro de 2026. **A correspondência central da janela LOC pendente ficou sustentada por alinhamento recíproco, permitindo consulta suplementar de variantes. Os flancos apresentam complexidade adicional. A investigação dos alinhamentos NPNT esclareceu que são alternativas imperfeitas, não cópias exatas equivalentes. Nenhuma janela foi promovida a safe harbor.**

## LOC: o que foi resolvido

Janela ROS NC_051811.1:17234956–17235956, 0-based half-open, correspondente ao intervalo UU NC_049228.1:17382280–17383266 (986 bp). O alinhamento original tem CIGAR `452M16I183M2D349M`, MAPQ60 e NM23. A decomposição conferida contra os FASTAs é de cinco substituições, 16 bases ROS sem par UU e duas bases UU sem par ROS.

O alinhamento da sequência UU de volta a ROS regressa exatamente às fronteiras da janela de 1 kb, MAPQ60, sem alternativa reportada, CIGAR `452M16D183M2I349M`. Os **984 pares de coordenadas são idênticos nos dois sentidos**, sem pares discordantes. Isto sustenta a correspondência desta janela central. Não demonstra igualdade entre as duas sequências nem genótipo de um animal experimental.

A regra operacional anterior NM≤10 impediu a consulta automática, apesar desta correspondência central. Foi mantida no resultado histórico, mas deixou de ser tratada como impedimento absoluto à obtenção de informação: fizemos uma consulta suplementar identificada como tal. Não foi baixado um filtro biológico nem alterado o scorer global.

## LOC: resultado suplementar de variantes

No intervalo UU correspondente, os [callsets públicos Dog10K](https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/) de 1.987 amostras devolveram:

| Callset | Registos | PASS | PASS com AF alternativa ≥1% |
|---|---:|---:|---:|
| SNP | 22 | 19 | 7 |
| nonSNP | 1 | 1 | 1 |

Todos os REF conferem com o FASTA UU. A leitura independente do texto VCF reproduz os totais, filtros e contagens de AF. Nenhuma SV sobreposta foi encontrada no catálogo SV consultado para esse intervalo; isso não prova ausência de SV.

Dos 23 registos, um SNP não tem projeção do intervalo REF contínua/disponível em ROS. Para os restantes, a tabela fornece apenas a projeção das coordenadas do REF. **Os alelos ALT não foram convertidos nem normalizados para ROS.** Mesmo um intervalo REF projetável não autoriza copiar automaticamente o alelo entre referências que diferem. As 16 bases ROS sem par UU continuam sem avaliação de SNPs por este VCF de referência UU; não são bases comprovadamente invariantes.

## LOC: os flancos não ficam validados

Também interrogámos a janela com 1 kb adicional de cada lado. Tanto no sentido ROS→UU como no inverso há alinhamentos divididos/suplementares. Os alinhamentos principais deixam 717 bases no início sem alinhar; no sentido inverso há ainda alternativas noutros contigs. Os SAMs completos foram guardados.

Assim, a reciprocidade exata da janela central **não se estende a toda a vizinhança de aproximadamente 3 kb**. A causa dos segmentos divididos pode envolver diferenças entre montagens e sequência repetitiva; não foi determinada nesta etapa. Este é um motivo para manter revisão de sequência local antes de qualquer desenho experimental. Não se inferiu uma alteração estrutural no dador.

## NPNT: natureza dos alinhamentos alternativos

Na janela NC_051836.1:27041661–27042661, os oito segmentos que falharam a proxy anterior abrangem, no conjunto, **NC_051836.1:27042361–27042586**. São quatro segmentos de 100 bp e quatro de 150 bp.

Foram registados oito alinhamentos primários exatos à origem e **31 alinhamentos alternativos reportados**. Estes últimos são parciais/imperfeitos, têm MAPQ máximo 0 e NM mínimo 0. Logo, o resultado anterior não demonstra cópias exatas igualmente boas. Demonstra que alguns segmentos têm semelhança com outros locais e falham a proxy conservadora de ausência de qualquer alternativa reportada.

A extensão indicada é o intervalo abrangido pelos segmentos testados, não a fronteira mínima de uma repetição. Retirar essas bases não foi demonstrado preservar o sinal T ou tornar a janela adequada. Não foi criada uma janela nova para evitar o aviso, nem realizada análise de especificidade de edição.

## Interpretação integrada

As 16 alternativas têm agora revisão local, incluindo consulta suplementar da única janela anteriormente não elegível, com limites de projeção explícitos. A janela LOC central deixou de ter correspondência por esclarecer, mas os seus flancos e as bases sem par permanecem pendências distintas. NPNT conserva o aviso de mapeabilidade, agora descrito sem confundir alternativas imperfeitas com duplicações exatas.

O resultado continua exploratório: o sinal T é escasso, vem de um único tumor e depende da definição celular. Ainda faltam QC/identidade completos do multiome, evidência replicada em T caninas relevantes, sequência individual e validação funcional após integração. O scoring das 461 regiões permanece v1.21.8; não foi reexecutado nesta etapa.

Ficheiros: `LOC_mapping_resolution.json`, `LOC_supplementary_lookup.json`, `LOC_supplementary_variants.tsv`, VCFs suplementares, `NPNT_failed_tile_alignments.tsv`, `NPNT_mapping_review.json`, `NPNT_alternative_summary.json`, `alternatives_integrated_v1222.tsv`, SAMs/logs de reciprocidade e flancos, `validation_summary.json`. As colunas históricas foram conservadas e a evidência nova acrescentada separadamente.
