from pathlib import Path
import json
N=Path(__file__).resolve().parent;A=N.parent/'ALL_CANDIDATES_20261005'
for p in [N/'whole_window_ROS_mappability_20261009/receipt.json',N/'ANO2_EpiC_basewise_review_20261009.json',N/'ANO2_original_regional_ATAC_20261009.json',N/'ANO2_conservation_identity_20261009_v3/receipt.json',N/'ANO2_indel_union_robustness_20261009.json',A/'integrated_evidence_20261008.json']:
 x=json.loads(p.read_text());print(p.name, list(x) if isinstance(x,dict) else 'list');print(str(x)[:4200] if 'integrated' not in p.name else str(next(r for r in x['results'] if '39708074' in str(r)))[:6500])
