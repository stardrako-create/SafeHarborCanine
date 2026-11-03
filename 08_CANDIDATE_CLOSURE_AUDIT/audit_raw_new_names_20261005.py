from pathlib import Path
import pysam,json
A=Path(__file__).resolve().parent;O=A/'W11_RAW_RECRUIT_20261005';names={x['name'] for x in json.loads((O/'comparison.json').read_text())['new_names_not_in_previous_local_selection']};rows=[]
with pysam.AlignmentFile(str(A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam'),'rb') as f:
 for x in f:
  if x.query_name in names:rows.append({'name':x.query_name,'mate':1 if x.is_read1 else 2,'ref':x.reference_name,'start0':x.reference_start,'cigar':x.cigarstring,'MAPQ':x.mapping_quality,'duplicate':x.is_duplicate,'qcfail':x.is_qcfail,'secondary':x.is_secondary,'supplementary':x.is_supplementary})
(O/'new_names_bam_audit.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows))