from pathlib import Path
import csv,json,shutil,hashlib,datetime
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'06_v1.21.1'
f=O/'cpg_islands_independent_calls.bed'
assert f.exists(), 'CpG build not finished'
with f.open() as inp:rows=list(csv.DictReader(inp,delimiter='\t'))
counts={}
for r in rows:
 t=r['type'];counts[t]=counts.get(t,0)+1
 assert t in ('CGI_GGF','CGI_TJ')
 minimum,gc,oe=(200,.5,.6) if t=='CGI_GGF' else (500,.55,.65)
 assert int(r['length'])==int(r['end'])-int(r['start'])>=minimum
 assert float(r['gc_frac'])>=gc and float(r['obs_exp'])>=oe
assert len({(r['#chrom'],r['start'],r['end'],r['type']) for r in rows})==len(rows)
(O/'cpg_validation.json').write_text(json.dumps({'call_counts':counts,'total':len(rows),'all_output_thresholds_pass':True,'duplicate_typed_coordinates':0},indent=2))
shutil.copy2(Path('/mnt/c/Users/Utilizador/OneDrive/Documents/New project/test_v1211.py'),O/'test_regression.py')
shutil.copy2(Path(__file__),O/'finish_run.py')
summary=json.loads((O/'summary.json').read_text())
report='''# Safe Harbor Canino — v1.21.1 corrigida

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
'''
for r in summary['comparison']:
 report+=f"| {r['genes']} | {r['previous_score']} | {r['corrected_score']} | {float(r['corrected_rrbs']):.2f}% | {float(r['rrbs_observed_fraction']):.2%} | {r['status']} |\n"
report+='''
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
'''
(O/'RESULTS.md').write_text(report,encoding='utf-8')
(O/'CURRENT.md').write_text('# v1.21.1 — scoring corrigido\n\nExecução concluída. Ler RESULTS.md. 1 candidato passa verificações registadas; 3 estão incompletos. Não são safe harbors validados. Usar code/ desta versão para repetir, não os scripts legados.\n',encoding='utf-8')
outputs={}
for name in ['candidates_scored_v1211.tsv','comparison.tsv','cpg_islands_independent_calls.bed','cpg_validation.json','RESULTS.md','test_regression.py']:
 p=O/name;outputs[name]=hashlib.sha256(p.read_bytes()).hexdigest()
(O/'output_hashes.json').write_text(json.dumps(outputs,indent=2))
print(json.dumps({'cpg':counts,'status':summary['counts'],'report':str(O/'RESULTS.md')},indent=2))
