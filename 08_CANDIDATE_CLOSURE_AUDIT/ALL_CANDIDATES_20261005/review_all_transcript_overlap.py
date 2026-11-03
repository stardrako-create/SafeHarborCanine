from pathlib import Path
import csv,json,collections,hashlib
O=Path(__file__).parent;P=O.parent.parent
rows=json.loads((O/'UU_projection/inputs.json').read_text())['windows']
rows=[dict(id=w['chrom']+'_'+w['start']+'_'+w['end'],chrom=w['chrom'],start0=w['start'],end0=w['end'],nearest_mRNA_TSS_distance=None,nearest_mRNA_id='',overlapping_transcripts=[],TAD_boundary_overlap=False) for w in rows]
bychrom=collections.defaultdict(list)
for r in rows:bychrom[r['chrom']].append(r)
features=collections.Counter()
with (P/'01_referencia/ROS_Cfam_1.0/genomic.gff').open() as f:
 for line in f:
  if line.startswith('#'):continue
  x=line.rstrip().split('\t')
  if len(x)!=9 or x[0] not in bychrom or not(x[2]=='mRNA' or x[2].endswith('RNA') or x[2]=='transcript'):continue
  features[x[2]]+=1;s=int(x[3])-1;e=int(x[4]);attrs=dict(y.split('=',1) for y in x[8].split(';') if '=' in y);identity=attrs.get('ID','')
  for r in bychrom[x[0]]:
   a=int(r['start0']);b=int(r['end0'])
   if s<b and e>a:r['overlapping_transcripts'].append(identity)
   if x[2]=='mRNA' and x[6] in ('+','-'):
    t=s if x[6]=='+' else e-1;d=max(t-b,a-(t+1),0)
    if r['nearest_mRNA_TSS_distance'] is None or d<r['nearest_mRNA_TSS_distance']:r['nearest_mRNA_TSS_distance']=d;r['nearest_mRNA_id']=identity
for line in (P/'04_tracks_processadas/ROS_Cfam_1.0/HiC/tad_boundaries.bed').open():
 if line.startswith('#') or not line.strip():continue
 x=line.split();s=int(x[1]);e=int(x[2])
 for r in bychrom.get(x[0],[]):
  if s<int(r['end0']) and e>int(r['start0']):r['TAD_boundary_overlap']=True
for r in rows:
 flags=[]
 if r['nearest_mRNA_TSS_distance'] is None:flags.append('mRNA_annotation_missing')
 elif r['nearest_mRNA_TSS_distance']<50000:flags.append('alternative_mRNA_TSS_under_50kb')
 if r['overlapping_transcripts']:flags.append('transcript_overlap')
 if r['TAD_boundary_overlap']:flags.append('TAD_boundary_overlap')
 r['additional_flags']=';'.join(flags);r['overlapping_transcripts']=';'.join(r['overlapping_transcripts']);r['pending']='TAD_risk_context;EpiC;repeats;conservation;variants;sequence_specificity;donor_genotype'
with (O/'transcript_and_boundary_review.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary=dict(reviewed=len(rows),additional_flags=dict(collections.Counter(r['additional_flags'] or 'none_in_this_layer' for r in rows)),feature_types_scanned=dict(features),status='additional_annotation_review_not_final_control_panel',source_gff_sha256=hashlib.sha256((P/'01_referencia/ROS_Cfam_1.0/genomic.gff').read_bytes()).hexdigest())
(O/'transcript_review_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
