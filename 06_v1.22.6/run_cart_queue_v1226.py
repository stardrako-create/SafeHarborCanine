from pathlib import Path
import subprocess,json,csv,time,shutil,traceback
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.6';B='/home/stardrako/miniforge3/envs/cart_rna/bin/';fastp='/home/stardrako/miniforge3/envs/atac/bin/fastp';manifest=list(csv.DictReader((P/'06_v1.22.5/raw_RNA_reprocessing_manifest.tsv').open(),delimiter='\t'));state={r['SRR']:{'status':'waiting_for_verified_input_and_index','sample':r['column']} for r in manifest};commands=[]
def save():
 tmp=O/'queue_status.tmp';tmp.write_text(json.dumps(dict(updated=time.time(),samples=state,commands=commands),indent=2));tmp.replace(O/'queue_status.json')
def execute(cmd,log):
 start=time.time()
 with log.open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 commands.append(dict(command=cmd,returncode=p.returncode,seconds=time.time()-start));save();assert p.returncode==0,(log,p.returncode)
def mark_complete(run):
 out=O/'samples'/run
 countfile=out/'STAR_ReadsPerGene.out.tab';lines=[x.split('\t') for x in countfile.read_text().splitlines()];assert len(lines)>4 and all(len(v)==4 and all(x.isdigit() for x in v[1:]) for v in lines)
 qc={}
 for l in (out/'STAR_Log.final.out').read_text().splitlines():
  if '|' in l:k,v=l.split('|',1);qc[k.strip()]=v.strip()
 gene_totals=[sum(int(v[i]) for v in lines[4:]) for i in [1,2,3]]
 state[run].update(status='completed_counts_pending_scientific_review',gene_rows=len(lines)-4,gene_count_totals_unstranded_forward_reverse=gene_totals,STAR_QC=qc)
# Recovery after a restart: `state` above always reinitializes every sample to
# 'waiting', but a sample whose STAR_ReadsPerGene.out.tab already exists and
# validates was genuinely finished before the process died - re-verify it and
# mark it complete instead of silently redoing fastp+STAR (each real run is
# ~10-45 min). A sample killed mid-STAR never wrote this file, so it correctly
# falls through and gets redone.
for r in manifest:
 if (O/'samples'/r['SRR']/'STAR_ReadsPerGene.out.tab').exists():
  try:mark_complete(r['SRR'])
  except Exception:pass
deadline=time.time()+7*86400;save()
try:
 while any(r['status']!='completed_counts_pending_scientific_review' for r in state.values()):
  if time.time()>deadline:raise RuntimeError('Seven-day bounded queue expired')
  ie=O/'reference/index_execution.json'
  if ie.exists() and json.loads(ie.read_text())['returncode']!=0:raise RuntimeError('STAR index generation failed; see index_execution.json')
  if not ie.exists():time.sleep(20);continue
  ready=None
  for r in manifest:
   if state[r['SRR']]['status']=='completed_counts_pending_scientific_review':continue
   paths=[O/'raw'/f"{r['SRR']}_{mate}.fastq.gz" for mate in [1,2]]
   if all(p.exists() and p.with_suffix(p.suffix+'.verified.json').exists() for p in paths):ready=(r,paths);break
  if not ready:
   ds=O/'download_status.json'
   if ds.exists() and any(v['status']=='failed' for v in json.loads(ds.read_text())['files'].values()):raise RuntimeError('An input download failed; queue paused without bypassing checksum')
   time.sleep(20);continue
  r,paths=ready;run=r['SRR'];out=O/'samples'/run;out.mkdir(parents=True,exist_ok=True);state[run]['status']='adapter_QC';save()
  assert shutil.disk_usage(O).free>30_000_000_000
  trimmed=[out/f'trimmed_{i}.fastq.gz' for i in [1,2]]
  cmd=[fastp,'--in1',str(paths[0]),'--in2',str(paths[1]),'--out1',str(trimmed[0]),'--out2',str(trimmed[1]),'--thread','4','--detect_adapter_for_pe','--disable_quality_filtering','--disable_length_filtering','--disable_trim_poly_g','--dont_eval_duplication','--json',str(out/'fastp.json'),'--html',str(out/'fastp.html')]
  execute(cmd,out/'fastp.log');fp=json.loads((out/'fastp.json').read_text());assert fp['summary']['before_filtering']['total_reads']>0
  state[run]['status']='STAR_alignment_and_gene_counts';save()
  # --outTmpDir must be on a native Linux filesystem: STAR opens a FIFO there
  # for its `zcat` readFilesCommand streaming, and DrvFs (WSL's /mnt/d mount
  # of the Windows D: drive) does not support FIFOs - STAR's own error names
  # this exact fix. Output files stay on D: via --outFileNamePrefix, which
  # doesn't need FIFO support. Per-run subdirectory: STAR refuses to reuse an
  # existing --outTmpDir, and the queue processes one sample at a time.
  star_tmp=f'/tmp/star_tmp_{run}'
  cmd=[B+'STAR','--runThreadN','6','--genomeDir',str(O/'reference/STAR_index'),'--readFilesIn',*[str(p) for p in trimmed],'--readFilesCommand','zcat','--outSAMtype','None','--quantMode','GeneCounts','--outFileNamePrefix',str(out/'STAR_'),'--outTmpDir',star_tmp]
  execute(cmd,out/'STAR_console.log');mark_complete(run);save();print(run,'counted',flush=True)
  # Keep raw and trimmed inputs for audit; no deletion is performed by this queue.
 # Only aggregate after all runs passed execution checks. Do not choose strandedness silently.
 allcounts={}
 for r in manifest:allcounts[r['column']]={v[0]:v[1:] for line in (O/'samples'/r['SRR']/'STAR_ReadsPerGene.out.tab').read_text().splitlines()[4:] if (v:=line.split('\t'))}
 names=[r['column'] for r in manifest];ids=list(allcounts[names[0]]);assert all(set(x)==set(ids) for x in allcounts.values())
 for col,label in enumerate(['unstranded','forward','reverse']):
  with (O/f'gene_counts_{label}.tsv').open('w') as f:
   f.write('gene_id\t'+'\t'.join(names)+'\n')
   for gene in ids:f.write(gene+'\t'+'\t'.join(allcounts[n][gene][col] for n in names)+'\n')
 (O/'raw_reanalysis_execution_complete.json').write_text(json.dumps(dict(status='execution_complete_requires_QC_strandedness_identity_and_gene_mapping_review',samples=9,genes=len(ids),limitations=['No differential-expression inference performed','Counts are STAR gene-assigned reads/fragments, not UMI-deduplicated molecules','All three count orientations preserved; orientation not inferred silently','Single-pass annotated ROS alignment, not replication of original CanFam3.1 pipeline']),indent=2));print('All9 counted; scientific review pending',flush=True)
except Exception as e:
 (O/'queue_failure.json').write_text(json.dumps(dict(error=str(e),traceback=traceback.format_exc(),time=time.time()),indent=2));raise
