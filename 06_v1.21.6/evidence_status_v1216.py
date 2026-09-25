from pathlib import Path
import csv,json
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.21.6'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
local=read(P/'06_v1.21.5/cohort_and_support.tsv');local={r['genes']:r for r in local if not r['genes'].endswith('_region')}
freq=read(O/'small_variant_frequencies.tsv');cons=read(O/'local_conservation_sensitivity.tsv');rows=read(P/'06_v1.21.5/local_windows_regulatory_rechecked.tsv')
matched=0
for r in rows:
 key=r['genes'];w=local[key];selected=r['start']==w['start'] and r['end']==w['end'];matched+=selected
 r['selected_1kb_previously']=str(selected)
 r['small_variant_review']='assessed_catalogue_variants_present' if selected else 'not_assessed'
 r['small_variant_PASS_SNP_records']=next(f['PASS_records'] for f in freq if f['genes']==key and '.SNPs.' in f['callset']) if selected else ''
 r['small_variant_PASS_nonSNP_records']=next(f['PASS_records'] for f in freq if f['genes']==key and '.nonSNPs.' in f['callset']) if selected else ''
 r['SV_review']=('large_catalogue_DEL_unresolved' if key.startswith('NPNT/') else 'no_overlapping_record_in_consulted_catalogue') if selected else 'not_assessed'
 r['conservation_status']='pilot_scale_sensitivity_measured' if selected else 'pilot_only'
 r['T_cell_accessibility_status']='not_validated'
 r['CAR_T_functional_safety_status']='not_validated'
 r['decision_status']='regulatory_conflict_requires_resolution' if int(r['epic_promoter_enhancer_bp_union'])>0 else 'insufficient_evidence_for_safe_harbor'
assert matched==3 and len(rows)==2091
with (O/'local_windows_evidence_status.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
selected=[r for r in rows if r['selected_1kb_previously']=='True']
(O/'evidence_status_validation.json').write_text(json.dumps(dict(rows=len(rows),selected_windows=matched,other_windows_variants_unassessed=sum(r['small_variant_review']=='not_assessed' for r in rows),no_validated_safe_harbor_claim=True,selected=[{k:r[k] for k in ['genes','start','end','epic_promoter_enhancer_bp_union','small_variant_PASS_SNP_records','small_variant_PASS_nonSNP_records','SV_review','decision_status']} for r in selected],rule='Evidence overlay, not a new ranking or automatic variant veto. Missing evidence remains unknown. One overlapping EpiC promoter/enhancer state is a regulatory conflict requiring review, not proof of in vivo activity.'),indent=2))
print(json.dumps(json.loads((O/'evidence_status_validation.json').read_text()),indent=2))
