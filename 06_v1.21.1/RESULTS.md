# Safe Harbor Canino — v1.21.1 corrigida

Data: 2026-09-12. Versão independente: conserva v120d e v1.21.0. O scoring foi executado sobre os 461 candidatos existentes e as tracks v1.20.0. Não é uma nova execução desde FASTQ, nem uma validação de safe harbors.

## Correções executadas

- RRBS: média apenas nas bases com evidência (coverage >=1 cão), aplicada igualmente a candidato e background de janelas. Zero medido é preservado; ausência de evidência é desconhecida. São médias de valores de bins ponderados por bases observadas, não uma média por CpG ou por reads.
- O builder RRBS desta versão foi também corrigido para escrever gaps na track média quando não existe evidência. Não foi preciso reconstruir as tracks legadas: o scorer usa a track de cobertura para as interpretar corretamente.
- Inputs de repetições, conservação, mappability, TAD e componentes do score ausentes produzem insufficient_evidence quando não existe outro veto. O BED de passing contém apenas passes_recorded_checks.
- Erros de leitura/intervalos inválidos propagam como erros, não são convertidos em zero. A amostragem de background termina com erro se não obtiver janelas observáveis suficientes.
- CpG: reconstrução genome-wide mantendo chamadas GGF e TJ separadas, com coordenadas e estatísticas próprias. Sobreposição não se transforma em falsa etiqueta de cumprimento simultâneo. Os limiares publicados são usados com o algoritmo de sliding/merge deste projeto; não alegamos equivalência exata a um programa publicado.
- Repetições: calculada a máscara soft da FASTA em todos os 461 candidatos. Difere até 16,92 pontos percentuais das 43 medições RepeatMasker existentes, por isso permanece diagnóstico separado, não substituição silenciosa.

## Resultado do scoring

457 excluídos, 3 com evidência insuficiente, 1 passa verificações registadas (ANO2/NTF3). Os 4 passam os vetos conhecidos, mas os 3 novos carecem de repeat_content e conservation. Não são 4 candidatos completamente avaliados.

| Candidato | Score anterior | Score corrigido | RRBS observado | Cobertura da janela | Estado |
|---|---:|---:|---:|---:|---|
| ANO2/NTF3 | 0.4309 | 0.5072 | 74.03% | 11.29% | passes_recorded_checks |
| NPNT/TBCK | 0.512 | 0.5015 | 81.49% | 8.38% | insufficient_evidence |
| LOC119876429/LOC119872513 | 0.4605 | 0.3988 | 84.71% | 6.16% | insufficient_evidence |
| LOC111090199/LOC100682550 | 0.3556 | 0.374 | 76.74% | 7.17% | insufficient_evidence |

A mediana do background RRBS corrigido é ~77,15%. ANO2/NTF3 tem ~74,03% nas bases observadas: o score relativo pode melhorar mesmo após corrigir uma média artificialmente baixa. Nenhum score é probabilidade de segurança, expressão ou sucesso. Não há novo veto arbitrário de metilação.

## Validação

8 testes de regressão passaram, incluindo escrita real de BigWig com gaps, zero observado, máscara do background, janela sem evidência, intervalo inválido e background sem amostras. O scoring completo terminou sem erro; validada a consistência dos estados de evidência. cpg_validation.json verifica todos os intervalos exportados contra os próprios limiares, tipos e duplicados.

## O que falta para candidatos biologicamente defensáveis

1. Completar repetições e conservação dos candidatos novos com referências e métodos comparáveis. O phyloP existente é um piloto com modelo neutro baseado numa região de 100 kb (ver V3_PROGRESS_NOTES.md); ter um valor não valida a qualidade desse modelo.
2. Comparar ATAC convencional e o agregado próprio, nos mesmos samples e intervalos, incluindo loci de Ehsan e picos reprodutíveis de referência. A normalização ou um percentil alto de médias de janelas não garante forte acessibilidade local.
3. Avaliar o tipo/estado celular pretendido: PBMC em conjunto não substitui evidência específica de células T caninas e do estado de ativação relevante. Medir consistência entre cães/cohorts e suporte local.
4. Definir o local concreto a avaliar dentro da região. O ponto de inserção e a vizinhança têm geometrias diferentes; não eliminar elementos regulatórios reais para forçar aprovação. Reconciliar o BED externo com assembly, versão, tipo e evidência. CpG é complementar e não valida atividade regulatória.
5. Separar segurança anotada, acessibilidade, metilação observada e incerteza. Avaliar ranks entre conjuntos fixos e resampling por cão/cohort; evitar percentagens inflacionadas por único sobrevivente.
6. Se nenhum dos 461 reunir as condições, expandir de forma predefinida o universo de busca, incluindo intervalos fora da restrição convergente de 50–75 kb. Isto exige nova geração e todos os mesmos controlos de risco; não relaxar limiares para recuperar nomes favoritos.
7. Validação experimental acordada com o grupo: integração correta e caracterizada, expressão estável sem depender de seleção contínua, ausência de alterações relevantes em genes locais/transcriptoma e função/viabilidade celular, consistência entre amostras independentes. Comparar com seleção aleatória e controlos mecanísticos; cassete/contexto celular fazem parte da conclusão.

Fontes: Ahmed et al. 2026 https://pmc.ncbi.nlm.nih.gov/articles/PMC12785581/ ; Aznauryan et al. 2022 https://pmc.ncbi.nlm.nih.gov/articles/PMC9017210/ ; Takai & Jones 2002 https://pmc.ncbi.nlm.nih.gov/articles/PMC122594/

## Ficheiros e execução

code/ contém a implementação corrigida usada nesta execução; manifest.json contém comando, hashes, runtime e inputs. candidates_scored_v1211.tsv guarda todos os candidatos; provisional_candidates.tsv conserva o grupo que passa os vetos conhecidos; candidates_passing_recorded_checks.bed exclui falta de evidência. comparison.tsv mostra as alterações. cpg_islands_independent_calls.bed é anotação suplementar, sem novo veto. run_scoring.py e test_regression.py permitem repetir o scoring/testes no ambiente atac existente.

Limites persistentes: tracks ATAC ainda usam o processamento anterior à recuperação de blocos corrompidos; conservação piloto; dados ausentes; anotação regulatória incompleta; média RRBS condicionada à cobertura e ao gate. Não foram desenhados gRNAs, enviados emails ou publicados resultados. Versões anteriores permanecem intactas.
