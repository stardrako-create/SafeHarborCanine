from pathlib import Path
import json,hashlib,collections
O=Path(__file__).parent
files=['UU_projection/inputs.json','UU_projection/projection_review.json','all_RNA_context.json','SV_all_windows/status.json','masked_specificity_20261008.json','pending_classification_20261006.json']
data={f:json.loads((O/f).read_text()) for f in files}
rows=data[files[0]]['windows'];mapping={r['id']:r for r in data[files[1]]['results']};rna={r['id']:r for r in data[files[2]]['results']};sv=data[files[3]];assert sv['state']=='completed'
seq={r['id']:r for r in data[files[4]]['results']};pending=data[files[5]]['windows'];out=[]
for w in rows:
 key=f"{w['chrom']}_{w['start']}_{w['end']}";hits=[h for h in sv['hits'] if h['query']==key]
 gaps=['T_cell_accessibility_not_validated','functional_safety_requires_experiments','complete_regulatory_isolation_not_demonstrated']
 if pending[key]['pending']:gaps.append('allele_reconciliation_or_spanning_deletion_context_pending')
 if mapping[key]['status']!='eligible_exact_window_for_variant_query':gaps.append('nonidentical_or_gapped_reference_projection')
 out.append(dict(id=key,genes=w['genes'],historical_decision=w.get('decision',w.get('review_disposition')),SV_catalogue_hits=hits,nearest_annotated_RNA_5prime=rna[key]['nearest'],exact_sequence_counts=seq[key]['counts'],repeat_bp=int(w['repeat_overlap_bp']),bulk_ATAC_percentile=float(w['atac_length_matched_percentile']),RRBS_observed_fraction=float(w['rrbs_observed_fraction']),RRBS_observed_mean=float(w['rrbs_mean_observed']),published_phyloP_max50bp_mean=w.get('published_phyloP_max50bp_mean'),queried_promoter_enhancer_bp=w['epic_promoter_enhancer_bp_union'],pending_projection=pending[key],remaining_evidence_gaps=gaps,disposition='not_advanced_unresolved_SV_context' if hits else 'computational_review_incomplete'))
assert len(out)==64
report=dict(results=out,counts=dict(collections.Counter(r['disposition'] for r in out)),approved_count=0,source_hashes={f:hashlib.sha256((O/f).read_bytes()).hexdigest() for f in files},limitations=['No new thresholds invented; integrated evidence ledger, not automatic biological approval','Historical exclusions remain applicable','Pending counts precede positional star-context review; positional overlap does not prove phase','No RNA boundary distance alone demonstrates regulatory safety'])
(O/'integrated_evidence_20261008.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report['counts']))