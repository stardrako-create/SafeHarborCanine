from pathlib import Path
import csv,json,gzip,collections,shutil
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.21.5';D=O/'external_sources'
with (P/'06_v1.21.4/all_local_windows.tsv').open() as f:local=list(csv.DictReader(f,delimiter='\t'))
chains=json.loads((P/'06_v1.21.4/chain_blocks.json').read_text());best={g:max(chains[g.split('/')[0]],key=lambda r:r['covered']) for g in {r['genes'] for r in local}}
annotations=[]
for file in sorted(D.glob('*_13_dense.bed.gz')):
 tissue=file.name.split('_')[0]
 with gzip.open(file,'rt') as f:
  for line in f:
   x=line.split()
   if int(x[3])==13:continue
   s,e=int(x[1]),int(x[2])
   for gene,ch in best.items():
    if x[0]!=ch['target']:continue
    assert ch['strand']=='+'
    for a,b,qa,qb in ch['blocks']:
     lo=max(s,qa);hi=min(e,qb)
     if lo<hi:annotations.append({'genes':gene,'tissue':tissue,'state':int(x[3]),'ros_start':a+lo-qa,'ros_end':a+hi-qa})
 print('Mapped regional annotation',tissue,flush=True)
def union(intervals):
 end=-1;total=0
 for s,e in sorted(intervals):total+=max(0,e-max(s,end));end=max(end,e)
 return total
def write(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write(O/'epic_nonquiescent_ROS_blocks.tsv',annotations)
for r in local:
 s,e=int(r['start']),int(r['end']);hits=[a for a in annotations if a['genes']==r['genes'] and a['ros_start']<e and a['ros_end']>s]
 reg=[a for a in hits if 1<=a['state']<=7]
 r['epic_nonquiescent_bp_union']=union([(max(s,a['ros_start']),min(e,a['ros_end'])) for a in hits])
 r['epic_promoter_enhancer_bp_union']=union([(max(s,a['ros_start']),min(e,a['ros_end'])) for a in reg])
 r['epic_regulatory_tissues']=';'.join(sorted({a['tissue'] for a in reg}))
 r['local_review_status']='regulatory_annotation_present' if reg else 'no_promoter_enhancer_in_queried_epic_tissues'
write(O/'local_windows_regulatory_rechecked.tsv',local)
summary=[]
for gene in sorted(best):
 rs=[r for r in local if r['genes']==gene and r['window_bp']=='1000'];top=max(rs,key=lambda r:float(r['atac_mean']))
 summary.append({'genes':gene,'windows_1kb':len(rs),'windows_with_promoter_enhancer':sum(r['epic_promoter_enhancer_bp_union']>0 for r in rs),'previous_top_1kb_flagged':top['epic_promoter_enhancer_bp_union']>0,'previous_top_1kb_regulatory_bp':top['epic_promoter_enhancer_bp_union'],'tissues':top['epic_regulatory_tissues']})
(O/'regulatory_recheck_summary.json').write_text(json.dumps(summary,indent=2));shutil.copy2(__file__,O/'regulatory_recheck_v1215.py');print(json.dumps(summary,indent=2))
