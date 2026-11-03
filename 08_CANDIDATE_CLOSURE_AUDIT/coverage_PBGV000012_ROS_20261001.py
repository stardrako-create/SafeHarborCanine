from pathlib import Path
import json,datetime,statistics,concurrent.futures,collections,pysam,fcntl
A=Path(__file__).resolve().parent;O=A/'SV_PBGV000012/ROS_structural_screen_20261001';B=A/'SV_PBGV000012/ROS_integration_20261001/selected_markdup.bam';R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fa';C='NC_051836.1'
def good(r):return not(r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail) and r.mapping_quality>=30
def work(start):
 end=start+10000
 with pysam.AlignmentFile(str(B),'rb') as f:depth=list(map(sum,zip(*f.count_coverage(C,start,end,quality_threshold=20,read_callback=good))))
 with pysam.FastaFile(str(R)) as f:seq=f.fetch(C,start,end).upper()
 acgt=sum(seq.count(x) for x in 'ACGT');gc=(seq.count('G')+seq.count('C'))/acgt if acgt else None
 zone='inside' if start>=6612614 and end<=35142119 else 'outside' if end<=6212614 or start>=35542119 else 'boundary_excluded'
 return {'start0':start,'end0':end,'zone':zone,'mean_depth':statistics.mean(depth),'median_depth':statistics.median(depth),'covered_bases':sum(d>0 for d in depth),'acgt_fraction':acgt/10000,'GC':gc}
if __name__=='__main__':
 lock=open(O/'coverage.lock','w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 (O/'coverage_status.json').write_text(json.dumps({'state':'running','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'workers':6}))
 with pysam.FastaFile(str(R)) as f:length=f.get_reference_length(C)
 with concurrent.futures.ProcessPoolExecutor(6) as pool:rows=list(pool.map(work,range(0,length-10000+1,100000)))
 selected=[r for r in rows if r['zone']!='boundary_excluded' and r['acgt_fraction']>=0.99]
 ins=[r['mean_depth'] for r in selected if r['zone']=='inside'];outs=[r['mean_depth'] for r in selected if r['zone']=='outside']
 bins=collections.defaultdict(lambda:collections.defaultdict(list))
 for r in selected:bins[int(r['GC']*20)][r['zone']].append(r['mean_depth'])
 gc=[]
 for k,v in sorted(bins.items()):
  if len(v['inside'])>=3 and len(v['outside'])>=3 and statistics.median(v['outside'])>0:gc.append({'GC_bin_low':k/20,'inside_n':len(v['inside']),'outside_n':len(v['outside']),'median_depth_ratio':statistics.median(v['inside'])/statistics.median(v['outside'])})
 report={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'windows':rows,'summary':{'inside_n':len(ins),'outside_n':len(outs),'inside_median_window_mean':statistics.median(ins),'outside_median_window_mean':statistics.median(outs),'median_depth_ratio':statistics.median(ins)/statistics.median(outs),'GC_strata':gc},'limits':['10kb windows every100kb; not exhaustive coverage','Approximate ROS flank boundaries with200kb exclusion margins','Within-chromosome controls; no genome-wide normalization','GC strata descriptive, no mappability matching or calibrated copy-number inference','Primary nonduplicate QCpass MAPQ>=30 and BQ>=20; repeats can distort depth']}
 (O/'coverage_review.json').write_text(json.dumps(report,indent=2));(O/'coverage_status.json').write_text(json.dumps({'state':'completed','utc':report['utc']}))
 hits=json.loads((O/'joint_hits.json').read_text())['hits'];scores=[]
 for name,hs in hits.items():
  g=max([h['AS'] for h in hs if not h['ref'].startswith('junction_')],default=None)
  j=[h for h in hs if h['ref'].startswith('junction_') and any(s<=980 and e>=1020 for s,e in h['blocks'])]
  scores.append({'name':name,'best_genomic_AS':g,'junction_spanning_hits_20bp_each_side':len(j),'best_spanning_junction_AS':max([h['AS'] for h in j],default=None)})
 (O/'joint_interpretation.json').write_text(json.dumps({'rows':scores,'criterion':'contiguous alignment block crossing synthetic offset1000 with>=20bp each side','scope':'four selected reads single-read competition; not paired likelihood'},indent=2))
 print(json.dumps({'coverage':report['summary'],'joint':scores}))