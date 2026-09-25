from pathlib import Path
import csv,json,hashlib,shutil,urllib.parse
from mem_local import call
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'08_CANDIDATE_CLOSURE_AUDIT';vr=json.loads((O/'variant_regulatory_review.json').read_text());ex=json.loads((O/'exact_sequence_review.json').read_text());assert vr['all_REF_and_haplotypes_verified'] and ex['status']=='completed'
windows=list(csv.DictReader((P/'07_FINAL_CANDIDATES_2026-09-19/shortlist_evidence.tsv').open(),delimiter='\t'));genes={};features=[]
with (P/'01_referencia/ROS_Cfam_1.0/genomic.gff').open() as f:
 for l in f:
  if l.startswith('#'):continue
  v=l.rstrip().split('\t')
  if len(v)!=9:continue
  a={k:urllib.parse.unquote(x) for item in v[8].split(';') if '=' in item for k,x in [item.split('=',1)]}
  if v[2] in ['gene','pseudogene']:genes[a['ID']]=dict(chrom=v[0],start=int(v[3])-1,end=int(v[4]),strand=v[6],symbol=a.get('Name',a['ID']),biotype=a.get('gene_biotype',''))
  elif v[2] not in ['exon','CDS','region']:features.append((v,a))
context=[]
for w in windows:
 s,e=int(w['start']),int(w['end']);local=[g for g in genes.values() if g['chrom']==w['chrom']];gap=lambda g:max(g['start']-e,s-g['end'],0);nearest=min(local,key=gap)
 coding=[]
 for v,a in features:
  parent=genes.get(a.get('Parent',''))
  if v[0]!=w['chrom'] or not parent or parent['biotype']!='protein_coding' or v[6] not in ['+','-']:continue
  pos=int(v[3])-1 if v[6]=='+' else int(v[4])-1;coding.append(dict(gene=parent['symbol'],transcript=a.get('ID'),pos0=pos,gap_bp=max(pos-e,s-pos-1,0)))
 first=min(coding,key=lambda x:x['gap_bp']);context.append(dict(window_id=w['window_id'],nearest_all_gene=nearest,all_gene_body_gap_bp=gap(nearest),nearest_coding_transcript_5prime=first,coding_5prime_ge50kb=first['gap_bp']>=50000,all_gene_overlap=any(gap(g)==0 and g['start']<e and g['end']>s for g in local)))
(O/'coding_vs_all_gene_context.json').write_text(json.dumps(context,indent=2))
table='| Janela | Tamanho | Segmentos únicos | Segmentos com outras ocorrências exatas | Máximo de loci |\n|---|---:|---:|---:|---:|\n'
for r in ex['summary']:table+=f"| {r['window_id']} | {r['length']} | {r['unique']}/{r['segments']} | {r['nonunique']} | {r['max_locus_count']} |\n"
ct='\n'.join(f"- {r['window_id']}: corpo génico mais próximo considerando todos os biotipos = {r['nearest_all_gene']['symbol']} ({r['all_gene_body_gap_bp']:,} bp); 5′ de transcrito codificante mais próximo = {r['nearest_coding_transcript_5prime']['gene']} ({r['nearest_coding_transcript_5prime']['gap_bp']:,} bp)." for r in context)
note=f'''# Auditoria adicional de fecho — em curso, 19 setembro 2026

**O painel anterior é provisório, não encerra o trabalho solicitado.** O utilizador pediu continuação e updates de cinco minutos. Este checkpoint documenta trabalho novo, não apenas reclassificação dos resultados anteriores.

## Variantes resolvidas entre referências

Confirmada igualdade integral das duas sequências de 1 kb entre ROS e UU, incluindo reverse complement para w11. Convertidos 46 registos do catálogo, incluindo FILTER não-PASS preservados; normalização com bcftools1.24 produziu48 alelos. Todos os REF conferem com ROS e o haplótipo alternativo antes/depois da normalização é idêntico. A normalização recuperou bases minúsculas do FASTA softmasked; a validação é insensível a maiúsculas, mantendo a sequência. A versão inicial do validador parou nessa diferença de caixa e foi corrigida antes de aceitar os resultados.

Ficheiros: variant_projection_provenance.tsv, projected_unnormalized.vcf, projected_normalized_ROS.vcf, normalized_variant_alleles_ROS.tsv, variant_masks_PASS_ROS.tsv. AF guardada na tabela é o máximo por registo fonte, **não a AF específica de cada alelo desdobrado**. Não copiar esse máximo para uma conclusão de frequência alélica individual. As máscaras de REF não representam toda a área de efeito de um indel. Não foram inferidos genótipos individuais.

## Especificidade de sequência em todo o assembly

Percorridos {ex['bases_scanned']:,} bp em {ex['contigs']} contigs. Todos os segmentos consecutivos de20,25,50,100 bp nas duas janelas foram procurados exatamente nas duas orientações. Teste do algoritmo com sobreposições e fronteiras de blocos concordou com enumeração direta num exemplo controlado.

{table}

Isto fecha uma lacuna real da proxy antiga de100–250 bp: w01 tem19 segmentos de20 bp não únicos; w11 tem25 de20 bp e6 de25 bp. Não significa que todos sejam guias viáveis, que todos tenham PAM ou que sejam off-targets reais. Ainda falta análise de semelhanças com mismatches/bulges no contexto da nuclease/guia. Nenhuma sequência de guia foi escolhida. Contigs alternativos podem aumentar o número de ocorrências; todos os contigs do FASTA estão incluídos e os primeiros10 locais de cada padrão não único foram guardados.

## Regulação: corrigir o âmbito de “gene mais próximo”

O ficheiro histórico canine_all_genes_stranded.bed é de genes codificantes; os lncRNAs estavam noutro BED. As distâncias anteriores devem ser designadas **distâncias a genes codificantes naquele BED**, não distância a qualquer gene. GFF completo e RNA filho de gene foram agora usados para distinguir corpos e5′ por biotipo.

{ct}

O5′ anotado de lncRNA LOC119872513 fica a32.785 bp de w01; o5′ de LOC102156669 a30.888 bp de w11. São extremidades anotadas, não TSS medidos experimentalmente. A ausência de overlap lncRNA não prova isolamento regulatório.

Na revisão [Ahmed2026, secção2](https://pmc.ncbi.nlm.nih.gov/articles/PMC12785581/), o critério referido de50 kb usa o5′ de genes codificantes. Portanto, a presença desses lncRNAs próximos **não constitui automaticamente falha desse critério**. A checklist local confirma que o código histórico extraiu genes codificantes, embora a palavra “ALL” fosse usada de forma ambígua. Deve separar-se esse critério de proximidade de RNA não codificante. Um [estudo primário de2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10836832/) aplicou, entre outros filtros, distância superior a100 kb de lncRNA; isso é uma regra mais restritiva distinta, não um limiar que se possa atribuir silenciosamente ao nosso pipeline. Sob essa regra de sensibilidade, ambas as janelas atuais falhariam. Não se mudou o filtro para preservar ou excluir candidatos nesta etapa.

## Trabalho pendente ativo

1. Integrar posições normalizadas das variantes com os segmentos não únicos e explicitar subintervalos limitados por esses dados, sem os promover a novos safe harbors.
2. Avaliar especificidade com mismatches/PAM quando a nuclease estiver definida; foi perguntado ao utilizador qual nuclease e se existem sequências/genótipos dos cães.
3. Rever alcance/proveniência do contexto regulatório e SV, corrigindo afirmações excessivas, sem interpretar ausência de dados como aprovação.

O fecho experimental continua dependente do material individual e de evidência nas células relevantes. A automação de updates permanece **ativa**; não foi pausada pela produção deste checkpoint.
'''
(O/'CURRENT.md').write_text(note,encoding='utf-8')
for f in ['closure_variants_regulation.py','exact_sequence_closure.py','document_closure_progress.py']:shutil.copy2(f,O/f)
notes=Path(r'D:\Biblioteca\Biblioteca\Labs\Vasco Barreto Lab\Notas');dest=notes/'Safe Harbor CAR-T - Auditoria adicional de fecho - 2026-09-19.md';dest.write_text(note,encoding='utf-8')
# Retain earlier report as history, but clearly withdraw its claim of computational closure.
old=P/'07_FINAL_CANDIDATES_2026-09-19/LEIA_PRIMEIRO.md';text=old.read_text(encoding='utf-8');banner='> ESTADO SUBSTITUÍDO: painel provisório em auditoria adicional. Distâncias do BED anterior referem-se a genes codificantes, não todos os biotipos. Ver ../08_CANDIDATE_CLOSURE_AUDIT/CURRENT.md.\n\n'
if not text.startswith('> ESTADO SUBSTITUÍDO'):shutil.copy2(old,O/'previous_panel_brief.md');old.write_text(banner+text,encoding='utf-8')
(P/'CURRENT.md').write_text('# Estado atual: auditoria adicional de candidatos em curso\n\n[Checkpoint](08_CANDIDATE_CLOSURE_AUDIT/CURRENT.md). Painel anterior provisório. Variantes projetadas/normalizadas; unicidade exata20–100bp examinada; contexto coding versus lncRNA esclarecido. Updates de5min ativos.\n',encoding='utf-8')
r=call('mempalace_add_drawer',{'wing':'Locus Canino','room':'CANDIDATE_CLOSURE_AUDIT_ACTIVE_20260919','content':note,'source_file':str(dest),'added_by':'Codex'});assert not r.get('error') and not r.get('result',{}).get('isError');(O/'mempalace_progress_receipt.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(context,indent=2))
