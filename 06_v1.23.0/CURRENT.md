# v1.23.0 — auditoria de protocolo e estrutura dos reads

19 setembro 2026. **A descrição do kit e a estrutura dos FASTQs disponíveis não estão reconciliadas.** Esta etapa delimita a dúvida com dados; não identifica qual kit foi efetivamente usado nem invalida automaticamente a contagem anterior.

## Fontes verificadas

Descarregados os XML dos nove experimentos ENA ligados aos runs do manifesto. Todos declaram RNA-Seq, layout PAIRED e SMART-Seq v4 3’ DE. Essa concordância documental pode refletir metadados copiados do mesmo depósito; não é confirmação independente do procedimento experimental. Exemplo: [SRX22428478](https://www.ebi.ac.uk/ena/browser/api/xml/SRX22428478). A descrição concorda com o [GEO GSM7887650](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7887650).

O manual original da Takara, versão 092618, páginas 8 e 21, descreve R2 com um índice inline de seis bases seguido de poly(T), destinado a desmultiplexagem. Fonte primária de autoria do fabricante, consultada numa [cópia do manual](https://device.report/m/a1cff4f8333dfe571d210a1f84e734b6cdb0f611dc09ffa2d1e02c4b1109f9f6.pdf); a ligação atual do fabricante não ficou acessível nesta sessão. Esta descrição não deve ser extrapolada a todos os protocolos SMART-Seq nem aplicada cegamente a dados já processados.

## Verificação dos ficheiros disponíveis

Inspecionados os primeiros 100.000 reads de cada um dos 18 FASTQs ENA: 1,800,000 reads, amostra inicial determinística, não amostragem aleatória do ficheiro inteiro. Comprimentos encontrados: 150 bases. Valores Phred observados: 30.

Em R2, a fração com prefixo exatamente igual a um dos 12 índices do manual foi 0.186–0.222%. A fração com pelo menos 16 T nas posições 7–26 foi 0.069–0.101%. A assinatura conjunta índice+poly(T) ocorreu em 0.000–0.001% dos reads R2. O limiar 16/20 é uma definição exploratória explícita, não uma especificação do fabricante nem teste formal do kit.

Não se observa a assinatura esperada como estrutura dominante no início destes reads. Isto é compatível com processamento prévio, uma preparação diferente ou uma descrição incompleta; os dados não distinguem essas hipóteses. Não prova qual delas ocorreu. As contagens GeneCounts e o mapeamento elevado são evidência operacional adicional, mas não identificam o protocolo nem confirmam por si o papel de cada mate.

Consultados também campos submitted_format/ftp/md5/bytes da ENA para os nove runs. Ausentes os campos de formato e ligação de ficheiros originalmente submetidos em todas as respostas: True. Campos vazios não demonstram ausência dos originais nem estabelecem proveniência das qualidades constantes. Respostas brutas preservadas; não foram descarregados novos ficheiros de sequenciação.

## Decisão e limites

Manter a reanálise v1.22.7 como contagem descritiva dos pares disponibilizados, com esta ressalva explícita. Não remover seis bases automaticamente, não descartar R2 nem interpretar índices de amostra como UMIs sem suporte nos dados. Não promover inferência confirmatória de orientação/protocolo com base apenas no equilíbrio forward/reverse.

Uma análise de sensibilidade limitada comparando R1, R2 e pares pode testar se ambos os mates fornecem alinhamento génico consistente; ainda não foi executada nesta etapa. Não resolveria identidade física nem qual foi o kit. Para fechar proveniência experimental seria necessária documentação original da biblioteca/desmultiplexagem, ficheiros originais ou esclarecimento do depositante. Não foi enviado contacto a ninguém.

As conclusões anteriores de concordância de rótulos e presença de contagens mantêm-se como observações do material analisado. Não atribuir ausência génica a silêncio absoluto. Não há nova validação de acessibilidade ou segurança dos candidatos safe harbor.

## Artefactos

protocol_structure_audit.json, nove XML SRX, nove respostas submitted.tsv, protocol_v1230.py e finish_protocol_v1230.py. Sem alterações aos FASTQs, às contagens, aos filtros ou ao ranking. A recontagem concluída permanece v1.22.7; esta versão é uma auditoria, não uma nova execução do pipeline.
