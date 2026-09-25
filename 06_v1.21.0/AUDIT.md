# Auditoria independente — v1.21.0, 2026-09-12

Execução do scoring sobre as tracks existentes v1.20.0, mantendo p55, vetos e pesos. Esta versão acrescenta classificação explícita de evidência ausente, verificações de assembly/intervalos, hashes e diagnósticos locais. Não reconstruiu tracks nem gerou novos dados de conservação.

- Candidatos: 461.
- Estados de evidência: {'excluded': 457, 'insufficient_evidence': 3, 'passes_recorded_checks': 1}.
- Diferenças com v120d nos campos comparados: 0 (replay_comparison.json).
- Linhas CpG com etiqueta dupla que não satisfazem TJ nas coordenadas reportadas: 5635.

## Resultados materiais

1. Os três novos sobreviventes da v120d têm repetições e conservação por avaliar. A classificação anterior promovia ausência de dados a ausência de risco. Aqui ficam como insufficient_evidence, não como candidatos reprovados biologicamente.
2. Etiqueta CGI_GGF+CGI_TJ significa no código apenas sobreposição entre duas chamadas, não que o intervalo exportado satisfaça ambos os critérios. Preservar os dois intervalos e os seus identificadores/estatísticas é necessário.
3. 55% de intervalos externos a sobrepor CpG não se compara diretamente com a fração de bases do genoma em CpG. É necessário um modelo nulo ao nível dos intervalos; sobreposição não prova atividade de promotor.
4. As distâncias de 50 kb usam extremidades da janela; as de 300 kb usam o centro. São geometrias diferentes. A regra de 50 kb não define o futuro ponto de inserção. Atribuições a TADs são aproximações entre boundaries, não validação funcional.
5. O snapshot do builder não foi aplicado às tracks: estas ainda antecedem a recuperação de blocos corrompidos. Os diagnósticos aqui apresentados refletem as tracks existentes.
6. O relatório de sensibilidade antigo refere-se à v120c, não aos novos quatro da v120d. Percentagens de configurações não são probabilidades de segurança ou sucesso.

## Limites

passes_recorded_checks não significa safe harbor validado nem validação dos inputs de conservação/risco. Máscaras de cobertura não medem por si suporte biológico suficiente entre cães. CpG continua anotação complementar, sem novo veto. Não foram enviados emails, alteradas versões anteriores ou publicadas alterações.

## Reprodução

manifest.json contém o comando, inputs, hashes e runtime. scoring.log contém a saída completa. raw_scoring.tsv conserva a saída original; candidates_scored_v121.tsv acrescenta evidence_missing e audit_status. Os ficheiros foram separados para permitir auditar a mudança de interpretação sem alterar scores silenciosamente.

## Achado adicional prioritário: média RRBS confunde ausência de dados com zero

O builder escreve zero na média quando não há evidência. O scorer calcula a média sobre a janela inteira e apenas verifica se existe alguma cobertura algures na janela. A média reportada é por isso diluída pela fração de bases sem evidência. A estatística condicionada a coverage>0 foi medida diretamente nas quatro janelas e está em rrbs_missing_as_zero_audit.tsv. Esta é uma comparação diagnóstica, não uma recalibração do score: corrigir exige aplicar a mesma máscara aos candidatos e ao background. Todos os scores de low_methylation e rankings associados permanecem provisórios.

| locus | média atual (%) | fração com evidência | média nas bases com evidência (%) |
|---|---:|---:|---:|
| LOC119876429/LOC119872513 | 5.22 | 6.161% | 84.71 |
| LOC111090199/LOC100682550 | 5.50 | 7.166% | 76.74 |
| ANO2/NTF3 | 8.36 | 11.293% | 74.03 |
| NPNT/TBCK | 6.83 | 8.378% | 81.49 |
