from pathlib import Path
import subprocess,json,hashlib,re,time,collections
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.6';R=O/'reference';R.mkdir(exist_ok=True);B=Path('/home/stardrako/miniforge3/envs/cart_rna/bin');ref=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';gff=ref.parent/'genomic.gff';gtf=R/'ROS_release106.gtf';commands=[]
raw=R/'gffread_raw.gtf'
if not raw.exists() and gtf.exists() and gtf.stat().st_size>0:gtf.replace(raw)
cmd=[str(B/'gffread'),str(gff),'-T','-o',str(raw)]
if not raw.exists():
 with (R/'gffread.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 commands.append(dict(command=cmd,returncode=p.returncode));assert p.returncode==0
else:commands.append(dict(command=cmd,reused_previous_successful_conversion=True))
# NCBI direct gene-parent exons (e.g. pseudogenes) are emitted without gene_id by
# gffread. Restore only where transcript_id exactly identifies an original gene.
gene_features={}
for line in gff.open():
 if line.startswith('#'):continue
 v=line.rstrip().split('\t')
 if len(v)==9 and v[2] in ['gene','pseudogene']:
  a=dict(x.split('=',1) for x in v[8].split(';') if '=' in x);gene_features[a['ID']]=a.get('Name',a.get('gene',''))
repaired=0
with gtf.open('w') as out:
 for line in raw.open():
  v=line.rstrip().split('\t')
  if len(v)==9:
   a=dict(re.findall(r'(\w+) "([^"]*)"',v[8]))
   if 'gene_id' not in a:
    tid=a.get('transcript_id');assert tid in gene_features,(tid,'missing original gene parent');assert a.get('gene_name',gene_features[tid])==gene_features[tid]
    v[8]+=f' gene_id "{tid}";';line='\t'.join(v)+'\n';repaired+=1
  out.write(line)
(R/'direct_gene_parent_repair.json').write_text(json.dumps(dict(features_repaired=repaired,rule='Add gene_id only when transcript_id exactly matches original GFF gene ID and gene_name agrees; no exon coordinates altered'),indent=2))
lengths={v[0]:int(v[1]) for l in (ref.with_suffix('.fna.fai')).read_text().splitlines() if (v:=l.split('\t'))};genes=set();transcripts={};exons=0
for l in gtf.open():
 if l.startswith('#'):continue
 v=l.rstrip().split('\t');assert len(v)==9
 if v[2]!='exon':continue
 assert v[0] in lengths and 1<=int(v[3])<=int(v[4])<=lengths[v[0]];a=dict(re.findall(r'(\w+) "([^"]*)"',v[8]));assert 'gene_id' in a and 'transcript_id' in a;genes.add(a['gene_id']);old=transcripts.setdefault(a['transcript_id'],a['gene_id']);assert old==a['gene_id'];exons+=1
hashes={}
for p in [ref,gff,gtf]:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 hashes[p.name]=h.hexdigest()
versions={x:subprocess.check_output([str(B/x),'--version'],stderr=subprocess.STDOUT,text=True).strip() for x in ['STAR','gffread']}
(R/'annotation_validation.json').write_text(json.dumps(dict(status='passed',assembly='GCF_014441545.1 ROS_Cfam_1.0',annotation='NCBI Annotation Release106, existing GFF header',genes_with_exons=len(genes),transcripts=len(transcripts),exon_features=exons,gene_ids=sorted(genes),sha256=hashes,versions=versions,commands=commands),indent=2));print('Annotation checked',len(genes),'genes',exons,'exons',flush=True)
index=R/'STAR_index';index.mkdir(exist_ok=True)
cmd=[str(B/'STAR'),'--runMode','genomeGenerate','--runThreadN','6','--genomeDir',str(index),'--genomeFastaFiles',str(ref),'--sjdbGTFfile',str(gtf),'--sjdbOverhang','149','--genomeSAsparseD','2','--limitGenomeGenerateRAM','28000000000','--outFileNamePrefix',str(R/'index_')]
(R/'index_command.json').write_text(json.dumps(cmd,indent=2));start=time.time()
with (R/'index_console.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
manifest=dict(returncode=p.returncode,seconds=time.time()-start,command=cmd,versions=versions,reference_hashes=hashes);(R/'index_execution.json').write_text(json.dumps(manifest,indent=2));assert p.returncode==0
assert all((index/x).stat().st_size>0 for x in ['Genome','SA','SAindex','genomeParameters.txt']);print('STAR index complete',flush=True)
