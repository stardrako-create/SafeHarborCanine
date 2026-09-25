from pathlib import Path
import json,csv,urllib.request,urllib.parse,gzip,hashlib,concurrent.futures,subprocess,os,collections
import pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.5';D=O/'external_sources';D.mkdir(exist_ok=True)
with (O/'cohort_and_support.tsv').open() as f:windows=[r for r in csv.DictReader(f,delimiter='\t') if not r['genes'].endswith('_region')]
chains=json.loads((P/'06_v1.21.4/chain_blocks.json').read_text());mapped=[]
for r in windows:
 gene=r['genes'].split('/')[0];c=max(chains[gene],key=lambda x:x['covered']);assert c['strand']=='+'
 for a,b,qa,qb in c['blocks']:
  s=max(a,int(r['start']));e=min(b,int(r['end']))
  if s<e:mapped.append({'genes':r['genes'],'chrom':c['target'],'start':qa+s-a,'end':qa+e-a,'source_start':s,'source_end':e})
def write(p,rows):
 if not rows:return
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write(O/'local_canfam3_blocks.tsv',mapped)
repo='https://api.github.com/repos/snu-cdrc/dog-reference-epigenome/git/trees/main?recursive=1'
tree=json.load(urllib.request.urlopen(repo));commit=tree['sha'];(D/'epic_tree.json').write_text(json.dumps(tree))
entries=[x for x in tree['tree'] if x['path'].startswith('Data/04_chrom_states') and x['path'].endswith('.bed')]
def fetch(entry):
 path=entry['path'];url='https://raw.githubusercontent.com/snu-cdrc/dog-reference-epigenome/'+commit+'/'+urllib.parse.quote(path,safe='/');name=Path(path).name;target=D/(name+'.gz')
 valid=False
 if target.exists():
  try:
   with gzip.open(target,'rb') as f:existing=f.read()
   valid=len(existing)==entry['size'] and hashlib.sha1(b'blob '+str(len(existing)).encode()+b'\0'+existing).hexdigest()==entry['sha']
  except (EOFError,OSError):pass
 plain=D/name
 if not valid and plain.exists():
  existing=plain.read_bytes()
  if len(existing)==entry['size'] and hashlib.sha1(b'blob '+str(len(existing)).encode()+b'\0'+existing).hexdigest()==entry['sha']:
   with gzip.open(target,'wb') as f:f.write(existing)
   valid=True
 if not valid:
  with urllib.request.urlopen(url,timeout=45) as inp,gzip.open(target,'wb') as out:
   while True:
    b=inp.read(1024*1024)
    if not b:break
    out.write(b)
 digest=hashlib.sha256();hits=[]
 with gzip.open(target,'rt') as f:
  for line in f:
   digest.update(line.encode());a=line.split()
   if len(a)<4 or a[0].startswith('#'):continue
   for m in mapped:
    if a[0]==m['chrom'] and int(a[1])<m['end'] and int(a[2])>m['start']:
     hits.append({'genes':m['genes'],'tissue_code':name.split('_')[0],'state':a[3],'overlap_bp':min(m['end'],int(a[2]))-max(m['start'],int(a[1])),'canfam3_chrom':m['chrom'],'annotation_start':int(a[1]),'annotation_end':int(a[2])})
 print('EpiC complete',name,flush=True)
 return {'url':url,'sha256_decompressed':digest.hexdigest(),'local':str(target)},hits
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(fetch,entries))
hits=[h for _,hs in results for h in hs];write(O/'epic_local_hits.tsv',hits)
counts=collections.Counter((h['genes'],h['tissue_code'],h['state']) for h in [])
for h in hits:counts[h['genes'],h['tissue_code'],h['state']]+=h['overlap_bp']
write(O/'epic_local_states.tsv',[{'genes':g,'tissue_code':t,'state':s,'observed_bp':n} for (g,t,s),n in sorted(counts.items())])
(D/'epic_sources.json').write_text(json.dumps({'commit':commit,'sources':[s for s,_ in results],'paper':'https://doi.org/10.1126/sciadv.ade3399','limitations':'11 tissues; not purified T cells; mapped source blocks only'},indent=2))
# Remap 1 kb sequences to the independently assembled Dog10K reference.
ref=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna';fa=pysam.FastaFile(str(ref))
with (O/'local_windows.fa').open('w') as out:
 for i,r in enumerate(windows):out.write(f'>{i}\n{fa.fetch(r["chrom"],int(r["start"]),int(r["end"]))}\n')
uu=Path('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');bin='/home/stardrako/miniforge3/envs/atac/bin/bwa'
with (O/'local_windows_uu.sam').open('w') as out,(O/'local_windows_uu.log').open('w') as err:subprocess.run([bin,'mem','-a','-t','4',str(uu),str(O/'local_windows.fa')],stdout=out,stderr=err,check=True)
vcf=pysam.VariantFile('/mnt/d/Jin2024_work/dog10k_sv/SV-genotype-v2.merge.agg_only.08032022.vcf.gz');matches=[];variants=[]
with pysam.AlignmentFile(str(O/'local_windows_uu.sam'),'r') as sam:
 rec=list(sam.fetch(until_eof=True));lengths=dict(zip(sam.references,sam.lengths))
 for i,r in enumerate(windows):
  al=[a for a in rec if a.query_name==str(i)];primary=[a for a in al if not a.is_secondary and not a.is_supplementary];assert len(primary)==1;a=primary[0]
  alternatives=len(al)>1 or a.has_tag('XA');nm=a.get_tag('NM') if a.has_tag('NM') else None
  eligible=not a.is_unmapped and not alternatives and a.mapping_quality>=30 and a.query_alignment_length==1000 and nm<=10
  # Exact reference contig length match; refuse ambiguous aliases.
  aliases=[name for name,c in vcf.header.contigs.items() if not a.is_unmapped and c.length==lengths[a.reference_name]]
  eligible=eligible and len(aliases)==1
  vname=aliases[0] if len(aliases)==1 else ''
  found=[]
  if eligible:
   for variant in vcf.fetch(vname,a.reference_start,a.reference_end):
    if variant.start<a.reference_end and variant.stop>a.reference_start:
     found.append({'genes':r['genes'],'vcf_chrom':vname,'start':variant.start,'end':variant.stop,'id':variant.id,'type':str(variant.info.get('SVTYPE','')),'AF':str(variant.info.get('AF','')),'filter':';'.join(variant.filter.keys())})
  matches.append({'genes':r['genes'],'uu_contig':a.reference_name,'vcf_contig':vname,'start':a.reference_start,'end':a.reference_end,'mapq':a.mapping_quality,'NM':nm,'query_aligned_bp':a.query_alignment_length,'alternatives':alternatives,'eligible_for_variant_lookup':eligible,'sv_records_overlapping':len(found) if eligible else None});variants+=found
write(O/'dog10k_local_mapping.tsv',matches)
if variants:write(O/'dog10k_local_SV.tsv',variants)
(O/'dog10k_lookup.json').write_text(json.dumps({'matches':matches,'SV_records':variants,'reference':'UU_Cfam_GSD_1.0; alias accepted only for exact unique contig-length match','source':'https://kiddlabshare.med.umich.edu/dog10K/Manta-SV_2022-03-28/SV-genotype-v2.merge.agg_only.08032022.vcf.gz','limitations':'SV callset only; no SNP/indel clearance; no claim of variant absence outside eligible aligned windows'},indent=2))
print(json.dumps({'states':[{'genes':g,'tissue':t,'state':s,'bp':n} for (g,t,s),n in sorted(counts.items())],'mapping':matches,'variants':variants},indent=2),flush=True)
