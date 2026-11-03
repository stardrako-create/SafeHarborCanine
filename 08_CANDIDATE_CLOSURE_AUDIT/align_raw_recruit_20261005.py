from pathlib import Path
import subprocess,json,datetime
A=Path(__file__).resolve().parent;O=A/'W11_RAW_RECRUIT_20261005';R=A/'SV_PBGV000010/alignment_chunked/ROS_joint_competition_20260926/genome_plus_junctions.fa'
cmd=['/home/stardrako/miniforge3/envs/atac/bin/bwa','mem','-a','-t','4',str(R),str(O/'selected_R1.fastq'),str(O/'selected_R2.fastq')];(O/'joint_command.json').write_text(json.dumps(cmd))
with (O/'joint_pairs.sam').open('w') as out,(O/'joint_pairs.log').open('w') as err:subprocess.run(cmd,stdout=out,stderr=err,check=True)
(O/'joint_status.json').write_text(json.dumps({'stage':'completed_needs_review','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}))