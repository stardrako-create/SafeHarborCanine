# w01 — controlos e decisões para preparar os testes

## O que já foi revisto
Os conjuntos existentes respondem a perguntas diferentes. Mantê-los separados evita chamar controlo experimental a um fundo computacional.

| Conjunto | Utilização adequada | Limite / decisão |
|---|---|---|
| Fundo ATAC v1.21.8 | Posicionar a janela relativamente ao genoma, usando a mesma largura | 3.000 janelas por largura. Observadas: 2.969 de 500bp, 2.971 de 1kb, 2.979 de 2kb. Comprimentos conferidos. Falta de sinal fica como não observado. Este fundo não está emparelhado por GC/mapeabilidade nem seleciona loci para edição. |
| Reamostragem v1.22.3 | Comparar sinal nos núcleos T com outros núcleos do mesmo conjunto, ajustando profundidade | 155 núcleos T consenso; 2.000 reamostragens, seed42. Mesmo tumor/dador; não são réplicas biológicas independentes. |
| CD3D/CD3E | Verificação da identidade/contexto T e das coordenadas de TSS | Coordenadas conferidas no recibo v1.21.7. Não escolher estes genes como controlos de integração. |
| Controlos WGS de w11 | Comparação descritiva de cobertura | Revisão GC/mapeabilidade concluída separadamente. Não transferir estes intervalos para o painel CAR-T. |

Fonte das contagens: 06_v1.21.8/local_background_windows.tsv; método em refresh_local_v1218.py. Fonte da reamostragem: 06_v1.22.3/depth_review_summary.json. Fonte TSS: 06_v1.21.7/control_TSS_validation.json.

## Matriz proposta para o ensaio de w01
Esta matriz define comparações e decisões; as condições operacionais e critérios numéricos dependem do sistema que o laboratório escolher.

| Comparação | Pergunta | O que precisa de ficar comparável |
|---|---|---|
| Células de origem sem intervenção | Qual é o estado basal? | Mesmo cão e estado de ativação/cultura |
| Controlo do procedimento de entrega | Que alterações resultam do procedimento? | Mesmo cão, manipulação e avaliação; componentes exatos a definir com a estratégia de entrega |
| w01 editado vs controlos do mesmo cão | A integração é íntegra, a expressão persiste e a função mantém-se? | Cassete, estado celular, momento de avaliação e critérios definidos previamente |
| Comparador de integração já caracterizado pelo laboratório, se existir | Como se comporta w01 relativamente a um sistema conhecido? | Mesma cassete e caracterização equivalente; nenhum locus positivo canino foi confirmado nos ficheiros revistos |

Avaliar estrutura e número de cópias, alterações locais, expressão e heterogeneidade, genes vizinhos, viabilidade e função. Separar eficiência de entrega/edição dos resultados obtidos nas células com integração confirmada. Usar cão como unidade biológica; clones e medições técnicas ficam identificados dentro de cada cão. A dimensão da experiência depende da variabilidade e do objetivo, ainda por definir.

## Controlos genómicos random e low-ATAC
São painéis com perguntas distintas. Um fundo aleatório permite contextualizar desempenho; um painel de baixa acessibilidade testa a relação com expressão. Selecionar por sinal baixo condiciona essa comparação e não produz um fundo aleatório.
Ainda não está fechado um painel de loci de controlo para edição: cada locus precisa de revisão própria de variantes, genes/regulação, sequência e viabilidade de desenho com a nuclease real. Não reutilizar automaticamente os 26 PASS ou os controlos antigos. Os ficheiros de fundo auditados acima não preenchem este requisito.

## Três entradas necessárias para fechar o desenho
1. Nuclease e variante exatas, mais estratégia de integração/entrega.
2. Ficheiro ou sequência da cassete CAR e identificação da versão pretendida.
3. Cães/amostras pretendidos e dados de sequência disponíveis para confirmar o locus.

O histórico previous_MOC.md, linhas próximas de 318 e 361, regista pedido dos plasmídeos reais ao Vasco. O dossier atual LEIA_PRIMEIRO.md mantém estas decisões pendentes. Não inferir a escolha a partir dos plasmídeos de Drosophila ou do estudo Cao. A pergunta ao utilizador continua pendente; nenhum contacto externo foi enviado.

## Próxima decisão
w01: preparar desenho quando chegarem as três entradas; conservar a janela e os ficheiros de variantes existentes. w11: continuar condicionado à resolução estrutural. Este documento fecha a distinção de funções dos controlos existentes e propõe a matriz de ensaio; a seleção concreta dos loci de controlo e o desenho de edição permanecem abertos.
