from pathlib import Path
import json,csv,hashlib,datetime,collections
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=C/'integrated_review_20260923'
def read(name):return json.loads((C/name).read_text())
split=read('split_conservation_20260923/status.json');assert split['state']=='completed' and not split['errors'],split['state']
O.mkdir(exist_ok=False)
rows=list(csv.DictReader((C/'transcript_and_boundary_review.tsv').open(),delimiter='\t'))
tad=read('TAD_multiscale_review.json');reg={r['id']:r for r in read('gapped_projection_20260923/CanFam3_partial_roundtrip_states.json')['summary']};repeats={r['id']:r for r in read('repeatmasker_context_20260923_v2/review.json')['results']};exact={r['window_id']:r['unique'] for r in read('exact_sequences/exact_sequence_review.json')['summary'] if r['length']==20};uu={r['id']:r for r in read('gapped_projection_20260923/UU_review.json')['results']}
full=read('variant_projection_20260923/review.json')['per_window'];part=read('partial_variant_projection_20260923_v3/review.json')['per_window'];oldcon={r['id']:r for r in read('conservation_20260923/alignment_support_review.json')['results']};newcon=collections.defaultdict(list)
for r in split['completed']:newcon[r['id']].append(r)
result=[]
for row in rows:
 key=row['id'];r=reg[key];u=uu[key];flags=[]
 if key in tad['risk_hit_ids']:flags.append('TAD_proxy_risk_gene')
 if key in tad['missing_assignment_ids']:flags.append('TAD_proxy_missing_assignment')
 if row['additional_flags']:flags.append(row['additional_flags'])
 if r['promoter_enhancer_tissues_in_projectable_bases']:flags.append('EpiC_promoter_enhancer_in_observed_bases')
 if r['unresolved_bp']:flags.append('regulatory_projection_incomplete')
 if u['additional_reported_alignments']:flags.append('UU_multiple_reported_alignments')
 if key not in full:flags.append('UU_whole_window_not_exact')
 if key in oldcon:scored=oldcon[key]['bp_with_non_dog_ACGT'] if oldcon[key]['status'].startswith('scores_available') else 0
 else:
  assert key in newcon;scored=sum(v['supported_scored_bp'] for v in newcon[key]);assert scored<=1000
 if scored<1000:flags.append('conservation_supported_coverage_incomplete')
 out={'id':key,'group':row['group'],'chrom':row['chrom'],'start0':int(row['start0']),'end0':int(row['end0']),'TSS_distance_bp':int(row['nearest_mRNA_TSS_distance']),'TAD_proxy_risk':key in tad['risk_hit_ids'],'TAD_assignment_missing':key in tad['missing_assignment_ids'],'regulatory_reciprocal_bp':r['reciprocal_unique_bp'],'regulatory_alert_tissues':r['promoter_enhancer_tissues_in_projectable_bases'],'repeat_bp':repeats[key]['masked_bp_union'],'unique_20mers_out_of_981':exact[key],'UU_exact_bp':u['exact_bp'],'UU_other_alignments':u['additional_reported_alignments'],'full_window_variant_records':full[key]['allele_records'] if key in full else None,'partial_only_variant_records':part.get(key),'conservation_supported_scored_bp':scored,'flags':flags,'status':'reviewed_computational_evidence_not_final_control_approval'}
 result.append(out)
assert len(result)==20 and len({r['id'] for r in result})==20
sources=['transcript_and_boundary_review.tsv','TAD_multiscale_review.json','gapped_projection_20260923/CanFam3_partial_roundtrip_states.json','repeatmasker_context_20260923_v2/review.json','exact_sequences/exact_sequence_review.json','gapped_projection_20260923/UU_review.json','variant_projection_20260923/review.json','partial_variant_projection_20260923_v3/review.json','conservation_20260923/alignment_support_review.json','split_conservation_20260923/status.json']
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':result,'sources_sha256':{f:hashlib.sha256((C/f).read_bytes()).hexdigest() for f in sources},'limits':['No new pass/fail thresholds or final control-panel selection','Full- and partial-window variant records kept separate; original counts preserved','Repeats and exact uniqueness reported as measurements, not new cutoffs','Conservation scores without non-dog ACGT support excluded from coverage','Donor genotypes, guide/nuclease/cassette and experimental validation outstanding']}
(O/'review.json').write_text(json.dumps(receipt,indent=2))
lines=['# Controlos: evidência integrada','', 'As 20 janelas foram revistas. Esta tabela junta as medições e mantém as lacunas visíveis; não aprova um painel final. Variantes parciais só descrevem trechos exatos. Conservação conta apenas posições com score e bases comparativas de outra espécie.','', '| Controlo | Repetidos bp/1000 | 20mers únicos/981 | Regulatório bp/1000 | Tecidos com alerta | Variantes: janela inteira / parcial | Conservação com suporte bp/1000 |','|---|---:|---:|---:|---|---|---:|']
for r in result:
 variant=str(r['full_window_variant_records'])+' / —' if r['full_window_variant_records'] is not None else '— / '+str(r['partial_only_variant_records']) if r['partial_only_variant_records'] is not None else 'não resolvido'
 lines.append(f"| {r['id']} | {r['repeat_bp']} | {r['unique_20mers_out_of_981']} | {r['regulatory_reciprocal_bp']} | {', '.join(r['regulatory_alert_tissues']) or 'nenhum nas bases avaliadas'} | {variant} | {r['conservation_supported_scored_bp']} |")
lines+=['','## Restrições por janela','']
for r in result:lines.append('- '+r['id']+': '+('; '.join(r['flags']) or 'Sem estes alertas nas camadas avaliadas; seleção final e validação experimental pendentes.'))
lines+=['','## O que continua por resolver','', 'Bases sem projeção e alelos que atravessam diferenças/gaps; regiões com múltiplos alinhamentos; origem faseada dos marcadores de deleção. Qualquer avaliação por janela inteira deve manter estas lacunas. w01 continua prioritário e w11 mantém a dúvida estrutural já documentada. A escolha da nuclease, cassete e cães dadores continua necessária para o desenho específico.','', 'Os ficheiros de origem e hashes estão em review.json. Outputs anteriores foram preservados.']
(O/'LEIA_PRIMEIRO.md').write_text('\n'.join(lines));print(json.dumps([{'id':r['id'],'supported_conservation_bp':r['conservation_supported_scored_bp']} for r in result],indent=2))
