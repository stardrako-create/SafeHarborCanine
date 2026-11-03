from pathlib import Path
import pysam,json,datetime,statistics
A=Path(__file__).resolve().parent;O=A/'SV_PBGV000012/ROS_structural_screen_20261001';O.mkdir(exist_ok=True)
B=A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam'
chrom='NC_051836.1';regions={'UU_left_flank_ROS_region':(35322119,35362119),'UU_right_flank_ROS_region':(6382614,6442614)}
def good(r):return not(r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail) and r.mapping_quality>=30
pairs={};sa=[];counts={}
with pysam.AlignmentFile(str(B),'rb') as f:
 for label,(start,end) in regions.items():
  n=0
  for r in f.fetch(chrom,start,end):
   if not good(r):continue
   n+=1
   pairs.setdefault(r.query_name,[]).append({'region':label,'read1':r.is_read1,'start0':r.reference_start,'end0':r.reference_end,'mapq':r.mapping_quality,'reverse':r.is_reverse,'mate_start0':r.next_reference_start,'cigar':r.cigarstring})
   if r.has_tag('SA'):
    for item in r.get_tag('SA').strip(';').split(';'):
     c,pos,strand,cigar,mq,nm=item.split(',');pos=int(pos)-1
     for other,(s,e) in regions.items():
      if other!=label and c==chrom and s<=pos<e and int(mq)>=30:sa.append({'name':r.query_name,'anchor_region':label,'anchor_mapq':r.mapping_quality,'SA_mapq':int(mq),'SA_start0':pos,'SA_cigar':cigar,'anchor_cigar':r.cigarstring})
  counts[label]=n
 links=[]
 for name,rs in pairs.items():
  for i,x in enumerate(rs):
   for y in rs[i+1:]:
    if x['region']!=y['region'] and x['read1']!=y['read1'] and x['mate_start0']==y['start0'] and y['mate_start0']==x['start0']:links.append({'name':name,'ends':[x,y]})
 cov=f.count_coverage(chrom,27041161,27042161,quality_threshold=20,read_callback=good)
 depths=list(map(sum,zip(*cov)))
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sample':'PBGV000012','reference':'ROS_Cfam_1.0','bam':str(B),'regions0':regions,'eligible_primary_records':counts,'reciprocal_mate_links_both_MAPQ30':links,'SA_links_anchor_and_SA_MAPQ30':sa,'w11_depth':{'start0':27041161,'end0':27042161,'mean':statistics.mean(depths),'median':statistics.median(depths),'covered_bases':sum(x>0 for x in depths),'base_quality_min':20,'mapping_quality_min':30},'limitations':['Approximate flank projection search regions, not exact reciprocal breakpoint liftover','Mate links and SA records are separate evidence types; do not sum them','No broad deletion coverage ratio computed here','Read-name groups are not confirmed independent molecules','No genotype confirmation or safe-harbor functional validation']}
(O/'review.json').write_text(json.dumps(r,indent=2));print(json.dumps({'paired_links':len(links),'SA_records':len(sa),'w11_depth':r['w11_depth']}))