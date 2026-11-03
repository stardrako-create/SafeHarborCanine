from pathlib import Path
import json,csv,pysam,hashlib
N=Path(__file__).resolve().parent;P=N.parent.parent
maps={r['id']:r for r in json.loads((N/'EpiC_control_projection.json').read_text())['results']};rows={r['id']:r for r in csv.DictReader((N/'transcript_and_boundary_review.tsv').open(),delimiter='\t')}
ref=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';fa=pysam.FastaFile(str(ref));results=[]
for label in ['1059','1266','1900']:
 key='bg1k_'+label;w=next(x for x in json.loads((N/'next_queue.json').read_text()) if x['id']==key);rid=f"{w['chrom']}_{w['start0']}_{w['end0']}";m=maps[rid];r=rows[rid];f=N/'canFam3_reference_identity_20261009'/(label+'.response.json');src=json.loads(f.read_text());assert m['status']=='eligible_unique_contiguous_base_mapping';assert src['chrom']==m['target'] and src['start']==m['start0'] and src['end']==m['end0']
 dna=src['dna'].upper();dna=dna if m['strand']=='+' else dna.translate(str.maketrans('ACGT','TGCA'))[::-1];ros=fa.fetch(r['chrom'],int(r['start0']),int(r['end0'])).upper();assert len(dna)==len(ros)==1000
 mismatch=[dict(offset0=i,ROS_position0=int(r['start0'])+i,ROS=a,CanFam3=b,CanFam3_position0=m['start0']+i if m['strand']=='+' else m['end0']-1-i) for i,(a,b) in enumerate(zip(ros,dna)) if a!=b]
 results.append(dict(id=key,strand=m['strand'],matched_bp=1000-len(mismatch),mismatches=mismatch,ROS_sequence_sha256=hashlib.sha256(ros.encode()).hexdigest(),CanFam3_oriented_sha256=hashlib.sha256(dna.encode()).hexdigest(),response_sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
out=N/'ROS_CanFam3_reference_identity_20261009.json';assert not out.exists();out.write_text(json.dumps(dict(results=results,limitation='Coordinate mapping is complete and round-trip reviewed; assembly sequence differences explicitly retained, not treated as individual variants or biological confirmation.'),indent=2));print(json.dumps(results))