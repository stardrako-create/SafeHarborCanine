from pathlib import Path
import json,pysam,hashlib,datetime
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=C/'reference_allele_reconciliation_20260923';d=json.loads((O/'review.json').read_text());n=0;covered=0
for row in d['records']:
 if row['status']!='allele_set_rebased_counts_preserved':continue
 f=C/'variant_catalogue_partial_query_20260923'/row['source_file']
 with pysam.VariantFile(str(f)) as vf:r=next(r for i,r in enumerate(vf) if i==row['source_record_index0'])
 originalcounts=[r.info['AN']-sum(r.info['AC'])]+list(r.info['AC']);back=[None]*len(originalcounts)
 for rosidx,sourceidx in enumerate(row['source_indices_in_ROS_order']):back[sourceidx]=row['ROS_counts'][rosidx]
 assert back==originalcounts and sum(back)==row['original_AN']
 assert row['ROS_alleles'][0]==row['ros_reference'];assert row['original_AF']==list(r.info['AF'])
 n+=1;covered+=len(row['pending_ids'])
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_records_reread':n,'pending_ALT_entries_in_reconciled_records':covered,'checks':['Inverse permutation reproduces entire source AC/AN count vector','Original source AF reread unchanged','ROS allele index0 equals verified reference'],'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [O/'review.json',O/'normalized_alleles.json',O/'rebased_normalized.vcf']}}
(O/'validation.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
