from pathlib import Path
import json,csv,pysam,hashlib
N=Path(__file__).resolve().parent;S=N/'hypothetical_TTTV23_alternate_PAM_sensitivity_20261009.json';x=json.loads(S.read_text());K='NC_051831.1_39708074_39709074';rr=[r for r in x['results'] if r['window']==K and not r['ungapped_hits'] and not r['single_bulge_hits']];assert len(rr)==2
fa=pysam.FastaFile(str(N.parent.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa'));rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];robust=json.loads((N/'ANO2_hypothetical_TTTV23_PAM_indel_robustness_20261009.json').read_text());assert robust['state']=='completed' and all(not r['failed'] for r in robust['checks']);out=[]
for i,r in enumerate(rr):
 s=r['start0'];seq=fa.fetch(r['chrom'],s,s+23).upper();seq=seq if r['strand']=='+' else rc(seq);assert seq==r['sequence'];pp=fa.fetch(r['chrom'],s-4,s).upper() if r['strand']=='+' else rc(fa.fetch(r['chrom'],s+23,s+27).upper());assert pp==r['PAM']
 out.append(dict(hypothesis_id='ANO2_TTTV23_'+str(i+1),chrom=r['chrom'],protospacer_start0=s,protospacer_end0=s+23,strand=r['strand'],protospacer_5to3=seq,PAM_5to3=pp,model='TTTV/CTTV/TCTV/TTCV; ungapped<=3 or onebasebulge<=2substitutions; ROSreference',off_locus_hits_in_defined_model=0,local5PASS_indel_combinations_with_intact_target_and_PAM=32,guide_approved=False,actual_nuclease_selected=False))
outfile=N/'ANO2_TTTV23_REFERENCE_MODEL_CANDIDATES_20261009.tsv';assert not outfile.exists()
with outfile.open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');writer.writeheader();writer.writerows(out)
receipt=dict(state='completed',source_sha256=hashlib.sha256(S.read_bytes()).hexdigest(),tsv_sha256=hashlib.sha256(outfile.read_bytes()).hexdigest(),results=out,limits=x['limits'],experimental_validation=False);(N/'ANO2_TTTV23_REFERENCE_MODEL_CANDIDATES_20261009.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(out))
