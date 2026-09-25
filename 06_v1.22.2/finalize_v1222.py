from pathlib import Path
import json,csv,shutil,hashlib
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.22.2'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
m=json.loads((O/'LOC_mapping_resolution.json').read_text());q=json.loads((O/'LOC_supplementary_lookup.json').read_text());n=json.loads((O/'NPNT_mapping_review.json').read_text());assert m['reciprocal_pair_intersection']==984 and m['forward_only_pair_count']==m['reciprocal_only_pair_count']==0
assert all(c['returncode']==0 for c in m['commands']);assert len(m['base_substitutions'])==5 and m['ROS_unpaired_bp']==16 and m['UU_unpaired_bp']==2
records=read(O/'LOC_supplementary_variants.tsv')
for s in q['summary']:
 lines=(O/(s['callset']+'.LOC_supplementary.vcf')).read_text().splitlines();rr=[l.split('\t') for l in lines if not l.startswith('#')];assert len(rr)==s['records'];assert sum(r[6]=='PASS' for r in rr)==s['PASS'];common=0
 for r in rr:
  info=dict(x.split('=',1) for x in r[7].split(';') if '=' in x);afs=[float(x) for x in info.get('AF','.').split(',') if x!='.'];common+=r[6]=='PASS' and bool(afs) and max(afs)>=.01
 assert common==s['PASS_AF_ge_0_01']
alternatives=[r for r in n['alignment_records'] if r['secondary'] or r['supplementary']];assert len(alternatives)==31
np_summary=dict(failed_tile_footprint=[n['source_footprint_start'],n['source_footprint_end']],alternative_alignment_records=len(alternatives),alternative_reference_contigs=sorted({r['reference'] for r in alternatives}),minimum_alternative_NM=min(r['NM'] for r in alternatives),maximum_alternative_MAPQ=max(r['mapq'] for r in alternatives))
(O/'NPNT_alternative_summary.json').write_text(json.dumps(np_summary,indent=2));rows=read(P/'06_v1.22.1/alternatives_integrated_v1221.tsv')
for r in rows:
 r['v1222_supplementary_lookup']='not_performed_in_this_patch';r['v1222_mapping_comment']='no_new_mapping_evidence_in_this_patch';r['supplementary_SNP_PASS']='';r['supplementary_nonSNP_PASS']=''
 if r['window_id']=='w00':
  r['v1222_supplementary_lookup']='UU_span_queried_after_reciprocal_core_confirmation; original_automated_gate_still_failed';r['v1222_mapping_comment']='984_reciprocal_pairs;16_ROS_unpaired_bp;2_UU_unpaired_bp;5_substitutions;flanks_split_and_alternatives';r['supplementary_SNP_PASS']=19;r['supplementary_nonSNP_PASS']=1
 if r['window_id']=='w13':r['v1222_mapping_comment']=f"8_failed_tiles;31_partial_imperfect_alternative_alignments;source_extent_{n['source_footprint_start']}_{n['source_footprint_end']};not_equally_good_exact_copies"
with (O/'alternatives_integrated_v1222.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
note=f'''# Safe Harbor CAR-T — v1.22.2: correspondência LOC e detalhe de mapeabilidade NPNT

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

Na janela NC_051836.1:27041661–27042661, os oito segmentos que falharam a proxy anterior abrangem, no conjunto, **NC_051836.1:{n['source_footprint_start']}–{n['source_footprint_end']}**. São quatro segmentos de 100 bp e quatro de 150 bp.

Foram registados oito alinhamentos primários exatos à origem e **31 alinhamentos alternativos reportados**. Estes últimos são parciais/imperfeitos, têm MAPQ máximo {np_summary['maximum_alternative_MAPQ']} e NM mínimo {np_summary['minimum_alternative_NM']}. Logo, o resultado anterior não demonstra cópias exatas igualmente boas. Demonstra que alguns segmentos têm semelhança com outros locais e falham a proxy conservadora de ausência de qualquer alternativa reportada.

A extensão indicada é o intervalo abrangido pelos segmentos testados, não a fronteira mínima de uma repetição. Retirar essas bases não foi demonstrado preservar o sinal T ou tornar a janela adequada. Não foi criada uma janela nova para evitar o aviso, nem realizada análise de especificidade de edição.

## Interpretação integrada

As 16 alternativas têm agora revisão local, incluindo consulta suplementar da única janela anteriormente não elegível, com limites de projeção explícitos. A janela LOC central deixou de ter correspondência por esclarecer, mas os seus flancos e as bases sem par permanecem pendências distintas. NPNT conserva o aviso de mapeabilidade, agora descrito sem confundir alternativas imperfeitas com duplicações exatas.

O resultado continua exploratório: o sinal T é escasso, vem de um único tumor e depende da definição celular. Ainda faltam QC/identidade completos do multiome, evidência replicada em T caninas relevantes, sequência individual e validação funcional após integração. O scoring das 461 regiões permanece v1.21.8; não foi reexecutado nesta etapa.

Ficheiros: `LOC_mapping_resolution.json`, `LOC_supplementary_lookup.json`, `LOC_supplementary_variants.tsv`, VCFs suplementares, `NPNT_failed_tile_alignments.tsv`, `NPNT_mapping_review.json`, `NPNT_alternative_summary.json`, `alternatives_integrated_v1222.tsv`, SAMs/logs de reciprocidade e flancos, `validation_summary.json`. As colunas históricas foram conservadas e a evidência nova acrescentada separadamente.
'''
(O/'CURRENT.md').write_text(note,encoding='utf-8');report=note+'\n\n---\n\n# Contexto histórico até v1.22.1\n\nA atualização acima resolve a correspondência central LOC e acrescenta a consulta suplementar; pendências históricas devem ser lidas à luz dessa atualização.\n\n'+(P/'06_v1.22.1/RELATORIO_ATUALIZADO_2026-09-14.md').read_text(encoding='utf-8');reportfile=O/'RELATORIO_ATUALIZADO_2026-09-15.md';reportfile.write_text(report,encoding='utf-8')
for name in ['resolve_mapping_v1222.py','loc_variants_v1222.py','finalize_v1222.py']:shutil.copy2(Path(name),O/name)
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');name='Safe Harbor CAR-T - v1.22.2 correspondencia LOC e NPNT - 2026-09-15';dest=notes/(name+'.md');dest.write_text(note,encoding='utf-8');shared=notes/'Safe Harbor CAR-T - Relatorio atualizado - 2026-09-15.md';shared.write_text(report,encoding='utf-8')
moc=notes/'MOC-SafeHarbor Canino CAR-T.md';text=moc.read_text(encoding='utf-8');anchor='# 00 · MOC — Safe Harbor Canino para CAR-T';assert anchor in text;shutil.copy2(moc,O/'previous_MOC.md');moc.write_text(text.replace(anchor,anchor+f'\n\n> [!important] Estado mais recente: v1.22.2 — 15 setembro\n> [[{name}]]; [[Safe Harbor CAR-T - Relatorio atualizado - 2026-09-15]]. Correspondência central LOC sustentada e variantes suplementares consultadas; flancos complexos. Alternativas NPNT parciais/imperfeitas. Sem safe harbor validado. Entradas seguintes são histórico.\n',1),encoding='utf-8')
shutil.copy2(P/'CURRENT.md',O/'previous_project_CURRENT.md');(P/'CURRENT.md').write_text('# Estado atual: v1.22.2 — 2026-09-15\n\n[Relatório atualizado](06_v1.22.2/RELATORIO_ATUALIZADO_2026-09-15.md) · [Tabela integrada](06_v1.22.2/alternatives_integrated_v1222.tsv).\n\nCorrespondência central LOC confirmada reciprocamente; consulta suplementar realizada. Flancos complexos; alternativas de alinhamento NPNT detalhadas. Scoring global v1.21.8; nenhuma validação funcional de safe harbor.\n',encoding='utf-8')
receipt=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'v1.22.2_Reciprocal_LOC_and_NPNT_Mapping','content':note,'source_file':str(dest),'added_by':'Codex'});assert not receipt.get('error') and not receipt.get('result',{}).get('isError');(O/'mempalace_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'validation_summary.json').write_text(json.dumps(dict(status='passed',reciprocal_pair_agreement_bp=984,NM_decomposition=dict(substitutions=5,ROS_unpaired=16,UU_unpaired=2),supplementary_variant_records=23,raw_VCF_count_filter_AF_checks='passed',NPNT_alternative_records=31,limits='Central reference correspondence only; flanks not fully collinear; no donor genotype or functional validation'),indent=2))
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='output_hashes.json'};(O/'output_hashes.json').write_text(json.dumps(hashes,indent=2));assert reportfile.read_text(encoding='utf-8')==shared.read_text(encoding='utf-8');print('v1.22.2 completed, checked and saved to project, Obsidian, MemPalace.');print(json.dumps(np_summary))
