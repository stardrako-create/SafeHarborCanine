from pathlib import Path
import csv,json,hashlib,shutil
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.5';s=json.loads((O/'source_audit_summary.json').read_text());runs=json.loads((O/'ENA_raw_read_inventory.json').read_text());assert s['ENA_runs']==9 and not s['ENA_errors'] and s['fastq_bytes']==63892132810
design=[]
for r in runs:
 assert len(r['records'])==1;x=r['records'][0];assert x['sample_alias']==r['GSM'] and x['experiment_accession']==r['SRX'];paths=x['fastq_ftp'].split(';');sizes=x['fastq_bytes'].split(';');assert len(paths)==len(sizes)==2
 design.append(dict(GSM=r['GSM'],column=r['column'],reported_donor_label=r['column'].split('_')[0],condition='untransduced' if r['column'].endswith('_T_cell') else 'B7H3' if '_B7H3_' in r['column'] else 'B7H3_CXCR2',SRX=r['SRX'],SRR=x['run_accession'],BioSample=x['sample_accession'],read1_url='https://'+paths[0],read2_url='https://'+paths[1],compressed_bytes=sum(map(int,sizes)),identity_status='consistent_accession_links; source_tissue_conflicts_preserved_in_v1224; not_genotype_verified'))
with (O/'raw_RNA_reprocessing_manifest.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(design[0]),delimiter='\t');w.writeheader();w.writerows(design)
assert len({r['SRR'] for r in design})==9
note='''# Safe Harbor CAR-T — v1.22.5: desenho emparelhado, genes ausentes e dados brutos

16 de setembro de 2026. **O artigo confirma três cães com condições emparelhadas. A pesquisa oficial de sinónimos não recuperou os 15 genes ausentes da matriz processada. Os nove runs de RNA-seq bruto foram localizados e associados às amostras.** Esta etapa esclarece a origem dos dados; não modifica os scores ou valida candidatos.

## Desenho confirmado

O texto e a legenda da figura 6 do [artigo original](https://doi.org/10.1007/s00262-024-03642-4) descrevem RNA-seq de três cães e comparação emparelhada de CAR-T B7-H3, CAR-T B7-H3/CXCR2 e T não transduzidas. Assim, o desenho é três indivíduos × três condições, não nove indivíduos. O método específico de RNA indica colheita após dez dias de expansão após transdução; a descrição geral das culturas admite dias 10–14. Esta distinção fica agora registada.

A correspondência de B/E/M aos rótulos de réplica do GEO foi verificada na v1.22.4. A confirmação do desenho no artigo sustenta a interpretação emparelhada ao nível do estudo. Não prova a identidade física de cada tubo. Os quatro campos source/tissue contraditórios continuam explícitos; não foram corrigidos por suposição ou escondidos.

## Auditoria de nomes génicos

Foi obtida a tabela corrente [NCBI Gene para Canis familiaris](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Canis_familiaris.gene_info.gz). Para os 40 genes de contexto foram pesquisados símbolo, símbolo oficial de nomenclatura, sinónimos e identificador LOC associado ao GeneID. Só seria aceite recuperação automática com correspondência única de gene e de linha da matriz.

**Nenhum dos 15 genes ausentes foi recuperado por um sinónimo inequívoco na matriz CPM.** Entre eles permanecem ANO2, NTF3, CD247 e CD8B. O resultado fecha a hipótese simples de recuperar esses valores apenas renomeando linhas com os sinónimos correntes consultados; não determina se a ausência veio de anotação histórica, filtragem ou baixa expressão. A tabela NCBI atual não é necessariamente a versão de anotação dos autores.

As matrizes anteriores mantêm esses valores ausentes, não zero. Não foi inferido silêncio génico nem ausência de função regulatória.

## Disponibilidade real de dados

O diretório público de suplementos de GSE247355 consultado contém apenas a matriz CPM já usada. O material suplementar ligado ao artigo é um PDF de figuras; não identificámos ali uma matriz de contagens brutas. Os métodos referem STAR/CanFam3.1, HTSeq e processamento RUVseq/DESeq2, mas não fornecem nessa descrição o ficheiro/release exato de anotação ou a regra que produziu exatamente as 14.385 linhas da matriz disponibilizada. Não se pode prometer reprodução exata a partir destes metadados.

Os nove SRX foram resolvidos no [ENA](https://www.ebi.ac.uk/ena/browser/home), cada um com um SRR e dois FASTQs. A associação GSM→SRX→SRR/BioSample foi verificada, com sample_alias correspondente ao GSM. Os ficheiros somam **63.892.132.810 bytes, aproximadamente 63,9 GB decimais comprimidos**. O manifesto tem condições, rótulos emparelhados, accessions, URLs e tamanhos. Estes aliases confirmam ligações de registos, não resolvem independentemente os campos biológicos contraditórios.

Os FASTQs não foram descarregados nem reprocessados nesta etapa. Reconstruir contagens a partir deles é uma nova análise computacional substancial, não uma conversão algébrica válida de CPM para contagens inteiras. A matriz CPM não foi utilizada em testes que exigem contagens brutas.

## Trabalho concretamente preparado

`raw_RNA_reprocessing_manifest.tsv` é a entrada rastreável para eventual reprocessamento dos nove runs. Será necessário fixar referência/anotação e versões, verificar integridade e QC, produzir contagens por gene e rever identidade amostral antes de usar um modelo emparelhado. Uma análise com anotação ROS atual deverá ser rotulada como reanálise, não reprodução idêntica do pipeline CanFam3.1 original.

A utilidade seria recuperar contexto de genes ausentes e comparar condições com contagens apropriadas. Continuaria a não medir diretamente acessibilidade das janelas nem neutralidade de integração. A revisão atual dos candidatos permanece exploratória e a necessidade de validação em T/CAR-T relevantes não foi removida.

## Ficheiros e verificação

`gene_alias_audit.tsv`, `source_audit_summary.json`, `ENA_raw_read_inventory.json`, nove respostas ENA TSV, `raw_RNA_reprocessing_manifest.tsv`, `paper_method_evidence.json`, texto XML do artigo, inventários dos diretórios públicos e snapshot NCBI Gene com hashes. As nove associações e pares FASTQ foram verificados; nenhuma recuperação de gene foi inventada. Scoring global permanece v1.21.8; contexto descritivo RNA permanece v1.22.4 com desenho agora confirmado no artigo.
'''
(O/'CURRENT.md').write_text(note,encoding='utf-8');report=note+'\n\n---\n\n# Histórico anterior até v1.22.4\n\nO desenho de três cães emparelhados e a auditoria de nomes/dados brutos foram atualizados acima.\n\n'+(P/'06_v1.22.4/RELATORIO_ATUALIZADO_2026-09-16.md').read_text(encoding='utf-8');(O/'RELATORIO_ATUALIZADO_2026-09-16.md').write_text(report,encoding='utf-8')
for name in ['source_audit_v1225.py','finish_v1225.py']:shutil.copy2(Path(name),O/name)
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.5 desenho genes e reads - 2026-09-16';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');shared=notes/'Safe Harbor CAR-T - Relatorio atualizado - 2026-09-16.md';shutil.copy2(shared,O/'previous_Obsidian_report.md');shared.write_text(report,encoding='utf-8');moc=notes/'MOC-SafeHarbor Canino CAR-T.md';m=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in m;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(m.replace(anchor,anchor+f'\n\n> [!important] Estado mais recente: v1.22.5\n> [[{name}]] — três cães emparelhados confirmados no artigo; genes ausentes não recuperados por sinónimos; nove runs brutos localizados (63,9 GB), ainda não reprocessados. Sem alteração da shortlist.\n',1),encoding='utf-8')
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.22.5\n\n[Relatório](06_v1.22.5/RELATORIO_ATUALIZADO_2026-09-16.md) · [Manifesto de reads RNA](06_v1.22.5/raw_RNA_reprocessing_manifest.tsv).\n\nDesenho de três cães emparelhados confirmado. Quinze genes ausentes não recuperados por sinónimos. Nove runs brutos inventariados, não reprocessados. Scoring e shortlist inalterados; sem safe harbor validado.\n',encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.5_RNA_Design_Alias_Raw_Read_Audit','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'output_hashes.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='output_hashes.json'},indent=2));assert shared.read_text(encoding='utf-8')==report;print('v1.22.5 source audit finalized and saved to project, Obsidian and MemPalace.')
