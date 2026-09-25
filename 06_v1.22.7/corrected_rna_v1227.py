"""Sequential corrected recount; immutable old outputs; fail-closed checkpoints."""
from pathlib import Path
import csv,json,hashlib,subprocess,time,shutil,fcntl,traceback,os
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
OLD=P/'06_v1.22.6'; O=P/'06_v1.22.7'; O.mkdir(exist_ok=True)
lock=open('/tmp/safeharbor_v1227.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
STAR='/home/stardrako/miniforge3/envs/cart_rna/bin/STAR'
rows=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'))
gtf=OLD/'audit_2026-09-18/ROS_release106.STAR_compatible.gtf'
fasta=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'
idx=O/'STAR_index';idx.mkdir(exist_ok=True)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,x):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2));tmp.replace(p)
state={'version':'1.22.7','started':time.time(),'status':'reference_validation','samples':{},'commands':[]}
def save():state['updated']=time.time();state['free_D_bytes']=shutil.disk_usage(O).free;write(O/'status.json',state)
def execute(cmd,log):
 start=time.time()
 with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 record=dict(command=cmd,returncode=r.returncode,started=start,finished=time.time());state['commands'].append(record);save()
 if r.returncode:raise RuntimeError(f'Failed command: {log}, code {r.returncode}')
def space():
 if shutil.disk_usage(O).free<35_000_000_000:raise RuntimeError('D free space below 35 GB; paused without deleting unvalidated inputs')
 if shutil.disk_usage('/mnt/c').free<15_000_000_000:raise RuntimeError('C free space below 15 GB; WSL temporary storage requires review')
try:
 save();space()
 annotation=json.loads((OLD/'reference/annotation_validation.json').read_text());repair=json.loads((OLD/'audit_2026-09-18/STAR_annotation_parser_impact.json').read_text())
 assert digest(gtf)==repair['repaired_sha256'];assert digest(fasta)==annotation['sha256'][fasta.name]
 refsig={'fasta_sha256':annotation['sha256'][fasta.name],'gtf_sha256':repair['repaired_sha256'],'STAR_version':subprocess.check_output([STAR,'--version'],text=True).strip()}
 command=[STAR,'--runMode','genomeGenerate','--runThreadN','6','--genomeDir',str(idx),'--genomeFastaFiles',str(fasta),'--sjdbGTFfile',str(gtf),'--sjdbOverhang','149','--genomeSAsparseD','2','--limitGenomeGenerateRAM','28000000000','--outFileNamePrefix',str(O/'index_')]
 refsig['command']=command
 receipt=O/'index_validated.json'
 if receipt.exists():assert json.loads(receipt.read_text())['signature']==refsig
 else:
  state['status']='building_corrected_index';save();execute(command,O/'index_console.log')
  ids=[l.split('\t')[0] for l in (idx/'geneInfo.tab').read_text().splitlines()[1:]]
  assert len(ids)==len(set(ids))==42309 and set(ids)==set(annotation['gene_ids'])
  assert 'finished successfully' in (O/'index_console.log').read_text()
  ih={n:digest(idx/n) for n in ['Genome','SA','SAindex','geneInfo.tab','genomeParameters.txt']}
  write(receipt,dict(signature=refsig,index_sha256=ih,genes=len(ids),validated=time.time()))
 indexreceipt=json.loads(receipt.read_text())
 for n,h in indexreceipt['index_sha256'].items():assert digest(idx/n)==h
 ids=[l.split('\t')[0] for l in (idx/'geneInfo.tab').read_text().splitlines()[1:]]
 allcounts={}
 for r in rows:
  run=r['SRR'];space();out=O/'samples'/run;out.mkdir(parents=True,exist_ok=True)
  inputs=[OLD/'samples'/run/f'trimmed_{m}.fastq.gz' for m in [1,2]]
  state['status']='verifying_trimmed_input';state['active_sample']=run;save()
  signature={'reference':refsig,'inputs':[{'path':str(p),'bytes':p.stat().st_size,'sha256':digest(p)} for p in inputs],'sample':r['column']}
  cmd=[STAR,'--runThreadN','6','--genomeDir',str(idx),'--readFilesIn',*[str(p) for p in inputs],'--readFilesCommand','zcat','--outSAMtype','None','--quantMode','GeneCounts','--outFileNamePrefix',str(out/'STAR_'),'--outTmpDir',f'/tmp/star_v1227_{run}_{os.getpid()}']
  # PID-specific scratch is operational, excluded from reproducibility signature.
  signature['alignment_parameters']=cmd[:-1]
  success=out/'validated.json'
  if success.exists():assert json.loads(success.read_text())['signature']==signature
  else:
   state['status']='aligning';save();execute(cmd,out/'console.log')
  assert 'ALL DONE!' in (out/'STAR_Log.out').read_text() and 'finished successfully' in (out/'console.log').read_text()
  qc={k.strip():v.strip() for l in (out/'STAR_Log.final.out').read_text().splitlines() if '|' in l for k,v in [l.split('|',1)]}
  data=[l.split('\t') for l in (out/'STAR_ReadsPerGene.out.tab').read_text().splitlines()]
  assert [v[0] for v in data[:4]]==['N_unmapped','N_multimapping','N_noFeature','N_ambiguous']
  assert [v[0] for v in data[4:]]==ids
  assert all(len(v)==4 and all(x.isdigit() for x in v[1:]) for v in data)
  npairs=int(qc['Number of input reads']);fp=json.loads((OLD/'samples'/run/'fastp.json').read_text())
  assert 2*npairs==fp['summary']['after_filtering']['total_reads']
  assert all(sum(int(v[i]) for v in data)==npairs for i in [1,2,3])
  hashes={n:digest(out/n) for n in ['STAR_ReadsPerGene.out.tab','STAR_Log.final.out','STAR_Log.out','console.log']}
  if success.exists():assert json.loads(success.read_text())['output_sha256']==hashes
  else:write(success,dict(signature=signature,output_sha256=hashes,QC=qc,validated=time.time()))
  allcounts[r['column']]={v[0]:list(map(int,v[1:])) for v in data[4:]}
  state['samples'][run]={'status':'validated','sample':r['column'],'unique_pct':qc['Uniquely mapped reads %'],'input_pairs':npairs};save()
  # Inputs retained while space is adequate. Never delete inputs before validation.
 names=[r['column'] for r in rows]
 for i,label in enumerate(['unstranded','forward','reverse']):
  dest=O/f'gene_counts_{label}.tsv';tmp=dest.with_suffix('.tmp')
  with tmp.open('w') as f:
   f.write('gene_id\t'+'\t'.join(names)+'\n')
   for g in ids:f.write(g+'\t'+'\t'.join(str(allcounts[n][g][i]) for n in names)+'\n')
  tmp.replace(dest)
 changes=[]
 for r in rows:
  prev={v[0]:list(map(int,v[1:])) for l in (OLD/'samples'/r['SRR']/'STAR_ReadsPerGene.out.tab').read_text().splitlines()[4:] if (v:=l.split('\t'))}
  new=allcounts[r['column']];common=set(prev)&set(new)
  changes.append(dict(sample=r['column'],added_genes=len(set(new)-set(prev)),changed_existing_genes_unstranded=sum(prev[g][0]!=new[g][0] for g in common),absolute_count_change_unstranded=sum(abs(prev[g][0]-new[g][0]) for g in common),added_gene_counts_unstranded=sum(new[g][0] for g in set(new)-set(prev))))
 write(O/'comparison_with_v1226.json',changes)
 state['status']='execution_validated_scientific_review_pending';save();write(O/'execution_complete.json',dict(samples=9,genes=len(ids),completed=time.time(),matrix_sha256={p.name:digest(p) for p in O.glob('gene_counts_*.tsv')},limitations=['Sample identity conflicts remain','Strandedness and inference require review','RNA expression does not validate safe harbor accessibility or safety']))
except Exception as e:
 state['status']='failed';state['error']=str(e);state['traceback']=traceback.format_exc();save();raise
