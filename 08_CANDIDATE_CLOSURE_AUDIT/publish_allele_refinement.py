from pathlib import Path
import json,shutil
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'08_CANDIDATE_CLOSURE_AUDIT'
r=json.loads((O/'allele_frequency_review.json').read_text());assert r['status']=='completed' and r['all_alleles_uniquely_matched']
note='''
## Frequências por alelo verificadas — 19 setembro, atualização seguinte

Concluída correspondência única dos 48 alelos normalizados aos ALT originais por equivalência do haplótipo local. A nova tabela `normalized_variant_alleles_with_AF_ROS.tsv` preserva o máximo por registo, acrescentando AF específica, índice do ALT original e identidade do alelo na referência UU. As tabelas anteriores permanecem intactas.

- w01: 27 alelos totais, 25 PASS, 10 alelos PASS com AF >=1%. O ALT T de w01_1 tem AF 0,302%, enquanto o máximo do registo é 13,5%; usar o máximo neste alelo teria alterado incorretamente a sua classe >=1%.
- w11: 21 alelos totais, 19 PASS, 9 alelos PASS com AF >=1%. O ALT TTGCTTC de w11_19 tem AF 0,05033%, versus máximo do registo 0,7801%; ambos abaixo de 1%.

Não confundir contagens de alelos com contagens de registos anteriores. Nenhuma variante foi retirada. As máscaras de REF de todas as variantes PASS e os resultados de unicidade exata permanecem iguais, porque não dependiam do limiar de AF. As AF são do catálogo, não dos cães que serão usados.

Revisitada a evidência já existente sobre o SV de w11: os quatro portadores catalogados têm 10–26 SNPs heterozigóticos de alta qualidade nos pequenos intervalos amostrados dentro da deleção alegada de 28,9 Mb. Isso questiona uma perda constitutiva simples nesses locais, mas a amostragem de SNPs não mede copy number de forma não enviesada, não verifica breakpoints e não exclui rearranjos complexos/mosaicismo. Não há nova validação por reads ou ensaio ortogonal nesta atualização. O SV mantém-se como limitação não resolvida; não foi convertido em falso positivo.

Próximas dependências: definir nuclease para especificidade com PAM/mismatches; localizar dados individuais se existirem; determinar se existe evidência bruta adicional adequada para o SV. A acessibilidade e função em CAR-T continuam a exigir dados relevantes. Sem worker pesado ativo nesta atualização; acompanhamento continua ativo.
'''
f=O/'CURRENT.md';t=f.read_text(encoding='utf-8').replace('Falta refinar a AF por alelo normalizado.','AF por alelo normalizado concluída; ver atualização seguinte.');f.write_text(t+note,encoding='utf-8')
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');dest=notes/'Safe Harbor CAR-T - Auditoria adicional de fecho - 2026-09-19.md';dest.write_text(f.read_text(encoding='utf-8'),encoding='utf-8')
receipt=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'ALLELE_AF_REFINED_20260919','content':note,'source_file':str(dest),'added_by':'Codex'})
assert not receipt.get('error') and not receipt.get('result',{}).get('isError')
(O/'mempalace_allele_AF_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
for name in ['refine_allele_frequencies.py','publish_allele_refinement.py']:shutil.copy2(name,O/name)
print('Updated audit checkpoint, Obsidian and MemPalace; originals preserved.')
