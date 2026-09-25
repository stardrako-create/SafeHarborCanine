"""Descriptive chr32 read-base coverage; never assigns SV truth automatically."""
from pathlib import Path
import csv,json,statistics,collections,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O=P/'08_CANDIDATE_CLOSURE_AUDIT';A=O/'SV_PBGV000010/alignment_chunked';CHROM='NC_049253.1'

def add_block(totals,start,end,size=100000):
 if end<=start:return
 for i in range(start//size,(end-1)//size+1):
  totals[i]+=max(0,min(end,(i+1)*size)-max(start,i*size))

def main():
 done=json.loads((A/'completed.json').read_text());assert done['status']=='alignment_completed_pending_biological_review'
 bam=pysam.AlignmentFile(done['bam'],'rb');assert bam.has_index()
 rows=list(csv.DictReader((O/'SV_chr32_100kb_reference_windows.tsv').open(),delimiter='\t'))
 modes=[(mq,dup) for mq in [20,30] for dup in [False,True]]
 sums={key:[0]*len(rows) for key in modes};reads=collections.Counter()
 for r in bam.fetch(CHROM):
  reads['all_records']+=1
  if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail:continue
  reads['primary_nonQCfail']+=1
  for mq,exclude_dup in modes:
   if r.mapping_quality<mq or (exclude_dup and r.is_duplicate):continue
   for start,end in r.get_blocks():add_block(sums[(mq,exclude_dup)],start,end)
 out=[];summary=[]
 for (mq,dup),values in sums.items():
  records=[]
  for row,total in zip(rows,values):
   rec=dict(row,MAPQ_min=mq,exclude_marked_duplicates=dup,aligned_read_bases=total,mean_read_base_depth=total/int(row['bases']))
   records.append(rec);out.append(rec)
  # Full-size, mostly canonical windows only; this is a stated descriptive filter.
  eligible=[r for r in records if int(r['bases'])==100000 and float(r['canonical_fraction'])>=0.99]
  groups={label:[r['mean_read_base_depth'] for r in eligible if r['relation_to_catalogued_span']==label] for label in ['inside','outside']}
  med={label:statistics.median(v) if v else None for label,v in groups.items()}
  summary.append(dict(MAPQ_min=mq,exclude_marked_duplicates=dup,eligible_window_counts={k:len(v) for k,v in groups.items()},median_depth=med,inside_outside_median_ratio=med['inside']/med['outside'] if med['outside'] else None))
 with (A/'coverage_windows_descriptive.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
 # Flags near boundaries are review leads, not deletion-supporting read counts.
 boundaries=[]
 for label,pos in [('left',4707938),('right',33613308)]:
  c=collections.Counter()
  for r in bam.fetch(CHROM,pos-10000,pos+10000):
   if r.is_unmapped or r.is_secondary or r.is_qcfail or r.mapping_quality<30:continue
   c['records_MAPQ30']+=1
   c['supplementary']+=r.is_supplementary;c['with_SA_tag']+=r.has_tag('SA');c['marked_duplicate']+=r.is_duplicate
   c['paired_not_proper']+=bool(r.is_paired and not r.is_proper_pair)
  boundaries.append(dict(label=label,catalogue_boundary0=pos,counts=dict(c)))
 result=dict(status='descriptive_metrics_complete_pending_review',coverage=summary,boundary_flag_counts=boundaries,record_counts=dict(reads),limitations=['Aligned read-base depth counts overlapping mates twice; not fragment depth','No base-quality filtering or GC/mappability correction in these metrics','Outside windows are not independently confirmed diploid controls','Duplicate flags derived from regional-retention BAM','SA/discordant flags do not establish breakpoint support or SV truth','Single catalogue carrier does not establish prospective CAR-T donor genotype'])
 (A/'coverage_and_boundary_descriptive.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

if __name__=='__main__':main()
