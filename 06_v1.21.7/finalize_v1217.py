from pathlib import Path
import json,csv,hashlib,shutil
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.21.7'
ref=json.loads((O/'ncbiRefSeq_CD3_region.json').read_text());targets=json.loads((O/'fragment_targets.json').read_text());control=[]
for g in ['CD3D','CD3E']:
 r=next(t for t in targets if t['name']==g+'_feature_control');mid=(r['start']+r['end'])//2;tx=[t for t in ref['ncbiRefSeq'] if t['name2']==g];assert tx and all((t['txStart'] if t['strand']=='+' else t['txEnd']-1)==mid for t in tx)
 control.append(dict(gene=g,TSS_position_0based=mid,feature_midpoint_matches_RefSeq_TSS=True,transcripts=[t['name'] for t in tx]))
(O/'control_TSS_validation.json').write_text(json.dumps(dict(source='https://api.genome.ucsc.edu/getData/track?genome=canFam6;track=ncbiRefSeq;chrom=chr5;start=15300000;end=15600000',downloadTime=ref['downloadTime'],controls=control),indent=2))
rows=list(csv.DictReader((P/'06_v1.21.6/local_windows_evidence_status.tsv').open(),delimiter='\t'))
for r in rows:
 if r['selected_1kb_previously']=='True':r['T_cell_accessibility_status']='exploratory_single_tumor_marker_masks_sparse_fragments_not_validated'
with (O/'local_windows_evidence_status.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
assert sum(r['T_cell_accessibility_status'].startswith('exploratory') for r in rows)==3
note='''# Checkpoint v1.21.7 — consulta local do multiome canino

2026-09-14. Trabalho técnico adicional; relatório de apresentação continua adiado. Este checkpoint acrescenta evidência exploratória, não valida os loci em células T ou CAR-T.

Montagem: as três janelas de 1 kb foram mapeadas de ROS_Cfam_1.0 para canFam6 por dois percursos: ROS→CanFam3→CanFam6 e inversão da chain CanFam6→ROS. Os conjuntos de coordenadas de destino concordam exatamente em todas as 3 000 bases, sem gaps nas janelas. Em canFam6: ANO2 chr27:7243778–7244778; LOC119876429 chr7:17088186–17089186; NPNT chr32:15200069–15201069 (0-based, half-open). Aliases chr/RefSeq obtidos do ficheiro UCSC canFam6.chromAlias.txt. Esta concordância não verifica o genótipo do doador. Dados e todas as alternativas de chain preservados em mapping_provenance.json.

Células: a matriz pública tem 5 969 núcleos; 5 849 passam os limites de contagens/features RNA e ATAC do script dos autores. Não reproduzimos WNN/cluster 4/8 nem afirmamos identidade celular definitiva. Os grupos abaixo são máscaras exploratórias de expressão, fixadas sem olhar para ATAC nos candidatos:
- 93 núcleos: pelo menos dois de CD3D/CD3E/CD3G detetados após QC.
- 58 núcleos: os mesmos critérios mais ≥1 de CD247/CD2/LCK e menos de dois marcadores mieloides (LYZ/LST1/CSF1R) e B (MS4A1/CD79A/CD79B).
- 133 núcleos: ≥1 CD3 e ≥1 marcador de suporte, com as mesmas exclusões de marcadores.
Dropout, doublets e identidade celular continuam por resolver. TRAC não foi acrescentado por correspondência especulativa de símbolo. Os grupos não são todos aninhados. As regras e os barcodes estão em cell_marker_masks.tsv; resumo em cell_mask_summary.json.

Fragmentos: o download integral (~1,92 GB comprimidos) estava lento e foi interrompido, preservando o ficheiro .partial. O ficheiro é BGZF ordenado por coordenadas. Fizemos pesquisa por blocos comprimidos via HTTP Range, usando a ordem dos contigs indicada no cabeçalho, e recuperámos intervalos que enquadram cada janela. As respostas Range, tamanhos, CRCs BGZF e ordenação local foram verificados; os blocos completos foram descodificados novamente por HTSlib e reproduzem os registos recuperados. Não há índice oficial local nem download global completo. A completude da consulta de inícios depende da ordenação do ficheiro; todas as janelas estão enquadradas por registos anteriores e posteriores. Contagens de sobreposição podem omitir fragmentos extremamente longos iniciados antes do intervalo recuperado, pelo que o resultado principal conta inícios dentro da janela.

Resultados: número de registos únicos de fragmentos com início dentro de cada janela de 1 kb, sem multiplicar pela coluna de suporte PCR:

| Janela | Todos após QC (5 849) | Grupo 93 | Grupo restrito 58 | Grupo 133 |
|---|---:|---:|---:|---:|
| ANO2/NTF3 | 25 | 0 | 0 | 1 |
| LOC119876429/LOC119872513 | 43 | 2 | 1 | 1 |
| NPNT/TBCK | 44 | 0 | 0 | 1 |

Os 43 registos LOC em todos os núcleos correspondem a 41 núcleos; NPNT 44 a 43; ANO2 25 a 25. Nos grupos de marcadores as contagens são de núcleos distintos. Não fizemos teste de enriquecimento que trate núcleos do mesmo tumor como replicados biológicos; não calculámos CPM de fragmentos sem denominador genómico. A soma de contagens da matriz de picos é uma métrica distinta e está apenas registada como contexto de biblioteca.

Controlos: as janelas de 1 kb centradas nas features CD3D/CD3E correspondem a TSSs RefSeq na mesma montagem, confirmação independente guardada em control_TSS_validation.json e ncbiRefSeq_CD3_region.json. No grupo restrito há somente dois fragmentos iniciados em cada controlo. Isso reforça a limitação de informação local. CD247/LYZ permanecem controlos associados à feature sem confirmação de TSS nesta execução. A confirmação posterior de CD3D/E está no JSON separado; os ficheiros brutos conservam a designação inicial prudente do controlo.

A matriz de picos não tinha uma feature sobre ANO2, embora a consulta direta recupere fragmentos nessa janela. LOC e NPNT têm cada um uma feature de pico sobreposta. Contagens da feature inteira não equivalem a contagens limitadas à janela; uma feature ausente não deve ser preenchida com zero de acessibilidade. Artefactos: candidate_peak_matrix_counts.tsv, local_T_marker_fragment_counts.tsv, regional_fragment_queries.json, fragment_quantification_validation.json e range_cache/.

Interpretação: este tumor de um único Doberman não fornece suporte convincente para acessibilidade forte de qualquer janela em células T. Os zeros não demonstram cromatina fechada. Um ou dois fragmentos em LOC não estabelecem vantagem robusta. A hipótese de priorizar LOC pelas anotações regulatórias/variantes anteriores continua apenas exploratória; não fica validada por estes dados. Também não elimina os conflitos EpiC de ANO2/NPNT ou a revisão SV pendente. Estado T atualizado apenas para as três janelas em local_windows_evidence_status.tsv; as outras 2 088 continuam sem avaliação T específica e nenhuma passa a safe harbor validado.

Fontes primárias: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244116 ; https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7807442 ; https://github.com/rln0005/OSA_snMultiomeSeq ; https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/ ; https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/liftOver/ ; https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/canFam6.chromAlias.txt .

Próximo trabalho que pode mudar a decisão: identidade celular com reanálise completa/labels reproduzíveis e QC de ATAC adequado; dados de T caninas saudáveis/ativadas com replicação biológica; avaliar janelas alternativas sem conflitos regulatórios preservando dados em falta; completar reconstrução global ATAC; resolver conservação piloto e SV NPNT. A validação funcional continua necessária. Não alterámos limites para fazer passar candidatos nem criámos relatório de apresentação.
'''
(O/'CURRENT.md').write_text(note,encoding='utf-8');shutil.copy2(P/'CURRENT.md',O/'previous_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.21.7\n\n[Checkpoint técnico](06_v1.21.7/CURRENT.md). Consulta local multiome: sinal escasso em grupos exploratórios de marcadores T; nenhum candidato validado. Relatório adiado.\n',encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - Checkpoint tecnico v1.21.7 - 2026-09-14';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_obsidian_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!important] Checkpoint v1.21.7\n> [[{name}]]: fragmentos locais do multiome medidos; sinal escasso nos grupos exploratórios T. Nenhum locus validado.\n',1),encoding='utf-8')
for f in Path('.').glob('*v1217.py'):shutil.copy2(f,O/f.name)
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.21.7_Technical_Cleanup','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'file_hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='file_hashes.json' and not p.name.endswith('.partial')},indent=2));print('Saved v1.21.7 project, Obsidian and MemPalace.')
