# w11 — avaliação estrutural após montagem local (05/10/2026)

**Decisão: condicionado. A deleção simples de 28,9 Mb não foi confirmada; a estrutura alternativa local também não está completamente reconstruída.**

As duas amostras apresentam sequência concordante na junção. SPAdes produziu contigs de154bp(PBGV000010) e175bp(PBGV000012); os154bp são idênticos após orientar ambos. AlinhamentosROS eUU recuperam77bp no flanco local e múltiplas localizações da porção repetitiva, comMAPQ0. Nenhuma montagem contém o motivo exato da cauda e o flanco oposto. Isto limita a montagem, não demonstra ausência da cauda.

A hipótese de inserção repetitiva local é compatível com duas fronteiras separadas14bp, cauda homopolimérica ancorada e pares locais. O teste exato20+20bp/BQ20 sustenta a junção da cauda em1nome da primeiraamostra e3da segunda. Existe também suporte de referência contínua nas duasamostras. Estes números não são estimativas de genótipo, nem prova de comprimento ou duplicação14bp.

A cobertura ampla não mostra a redução esperada de uma deleção heterozigótica simples: primeiraamostra razões normalizadas próximas1; segunda razão descritiva1,013608. Scores da junção, genótipo catalogado e modelosBREAKPOINT/COVERAGE não substituem demonstração de perda de cópias ou ligação ao locus distante. A projeção direita e a origem da repetição continuam ambíguas.

## Limite do trabalho concluído
Montagens usam janelas locais e mates presentes nos BAMs retidos. Não constituem recrutamento exaustivo dos SRA originais. As repetições podem excluir leituras relevantes do subconjunto e permitir caminhos artificiais; o percursok31 foi rejeitado como montagem sustentada. Reanálises e assemblers das mesmas leituras não são evidência independente.

## Trabalho computacional ainda justificável
Uma tentativa delimitada de recrutamento a partir dos dados brutos, com sementes específicas da junção/flancos e expansão por pares, pode recuperar leituras excluídas pelo filtro cromossómico. Tem de preservar origem por amostra, quantificar a seleção e controlar recrutamento de repetições. Não continuar apenas variando scores ou somando evidência repetida. Depois desta tentativa, reavaliar se existe sequência que ligue fisicamente as extremidades; se não, documentar resolução insuficiente dos dados.

## Evidência que resolveria a estrutura
Sequência que atravesse toda a alteração e chegue a flancos únicos de ambos os lados, ou evidência independente equivalente que distinga inserção local de rearranjo distante. Ainda seria necessário validar função e segurança do locus nas células relevantes: resolver a estrutura não valida biologicamente um safe harbor. O cão experimental não é requisito para continuar a presente investigação computacional.

## Resultados locais
- W11_SPADES_REVIEW_20261004/alignment_review.json
- W11_SPADES_PBGV10_REVIEW_20261005/alignment_review.json
- W11_SPADES_PBGV10_REVIEW_20261005/cross_sample_junction.json
- W11_TAIL_REFERENCE_20261002/review.json
- W11_PATH_AUDIT_20261002/review.json

Inputs, contagens e outputs anteriores preservados; sem publicação automática.

### 05/10/2026 05:06 UTC — montagem ampliada idêntica; acompanhamento automático encerrado
Montagemampliada concluiu05:04:46UTC. comparison.json: conjuntos das29sequências decontigs exatamente iguais aosanteriores(0exclusivos de cada montagem); cabeçalhos/coberturas ehashFASTA diferem. Junção continua em175bp, semmotivoexato da cauda. Os2paresnovos não prolongam a reconstrução. Não repetir alinhamento genómico de sequênciasidênticas; resultadosROS/UU anteriores continuam aplicáveis àssequências.
Encerrada a sequência delimitada de análise local, montagemdedicada em2amostras e recrutamentointegral porsementes na segunda. Não é alegação de esgotamento de todos osalgoritmos/dados possíveis. Não há execução útil definida restante nesta estratégia com osinputs atuais; novas iterações automáticas tenderiam a repetir evidência. Automação resolver-w11-acompanhamento-hor-rio eliminada conforme instrução de encerrar quando deixar de haver trabalho útil executável. Nenhumprocessoativo.
Conclusão: deleção simples28,9Mb não confirmada; alteração repetitiva local compatível com extremidades/cauda/pares, sem alelocompleto. w11 permanece condicionado, não validado biologicamente. Para resolver: ligação independente entreflancos únicos que atravesse a alteração(comprimento/sequência completos), ou novosdados equivalentes; não exigir cãoexperimental para análise de amostrasarquivadas. Contagens/inputs/outputs preservados, sem publicação. CURRENT/Obsidian eavaliação atualizados; MemPalace pendente emfila local.
