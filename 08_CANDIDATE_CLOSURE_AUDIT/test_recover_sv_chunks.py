import importlib.util, tempfile, gzip, json, random, subprocess, time
from pathlib import Path
src=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/08_CANDIDATE_CLOSURE_AUDIT/recover_sv_chunks.py')
s=importlib.util.spec_from_file_location('worker',src);w=importlib.util.module_from_spec(s);s.loader.exec_module(w)
t=Path(tempfile.mkdtemp(prefix='sv_chunk_test_'));w.D=t;w.A=t/'alignment_chunked';w.REF=t/'ref.fa';w.SIZE=1
r=random.Random(45);seq=''.join(r.choice('ACGT') for _ in range(4000));w.REF.write_text('>'+w.CHROM+'\n'+seq+'\n')
subprocess.run([str(w.B/'bwa'),'index',str(w.REF)],check=True,capture_output=True)
rev=lambda x:x.translate(str.maketrans('ACGT','TGCA'))[::-1]
for mate in (1,2):
 f=t/f'SRR15734832_{mate}.fastq.gz'
 with gzip.open(f,'wt') as h:
  for i in range(2):
   start=500+i*500;read=seq[start:start+150] if mate==1 else rev(seq[start+200:start+350])
   h.write(f'@pair{i}/{mate}\n{read}\n+\n'+150*'I'+'\n')
 (t/(f.name+'.validated.json')).write_text(json.dumps({'validated':True,'bytes':f.stat().st_size}))
(t/'fastq_pair_validation.json').write_text(json.dumps({'pairs':2,'gzip_streams_read_to_EOF':True,'status':'validated_FASTQ_structure_and_mate_identity'}))
original_sleep=time.sleep;w.time.sleep=lambda n:original_sleep(.05)
w.main()
assert json.loads((w.A/'completed.json').read_text())['chunks']==2
for i in range(2):
 c=w.A/f'chunk_{i:03d}';v=json.loads((c/'validated.json').read_text());assert v['primary']==2 and v['read1']==v['read2']==1
 assert not (c/'replay.sam').exists() and not (c/'mate1.fastq').exists()
assert all((t/f'SRR15734832_{i}.fastq.gz').exists() for i in (1,2))
count=subprocess.check_output([str(w.B/'samtools'),'view','-c','-F','2304',str(w.A/'selected_markdup.bam')],text=True).strip();assert count=='4'
print(json.dumps({'status':'passed','test_dir':str(t),'chunks':2,'primary_records':4,'originals_preserved':True,'completed_chunk_temporaries_removed':True}))
