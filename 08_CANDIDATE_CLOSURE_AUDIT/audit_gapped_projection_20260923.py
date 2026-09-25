from pathlib import Path
import pysam,csv,json,gzip,collections,datetime
A=Path(__file__).resolve().parent;C=A/'CONTROL_LOCUS_REVIEW_20260922';O=C/'gapped_projection_20260923';O.mkdir(exist_ok=False)
rows={r['id']:r for r in csv.DictReader((C/'transcript_and_boundary_review.tsv').open(),delimiter='\t')}
R=A.parent/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'
UU='/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa';rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
results=[];allmaps={}
with pysam.AlignmentFile(str(C/'UU_projection/full_reference.sam')) as sam,pysam.FastaFile(str(R)) as ref,pysam.FastaFile(UU) as uu:
 grouped=collections.defaultdict(list)
 for r in sam:grouped[r.query_name].append(r)
 for key,w in rows.items():
  prim=[r for r in grouped[key] if not r.is_unmapped and not r.is_secondary and not r.is_supplementary]
  assert len(prim)==1;r=prim[0];seq=ref.fetch(w['chrom'],int(w['start0']),int(w['end0'])).upper();oriented=rc(seq) if r.is_reverse else seq
  assert r.query_sequence.upper()==oriented
  m=[None]*1000;mismatch=0;exact=0;refonly=0;queryonly=0
  for q,p in r.get_aligned_pairs():
   if q is None:refonly+=1;continue
   if p is None:queryonly+=1;continue
   b=uu.fetch(r.reference_name,p,p+1).upper();same=b==oriented[q] and b in 'ACGT';i=999-q if r.is_reverse else q
   m[i]={'uu_pos0':p,'exact':same};exact+=same;mismatch+=not same
  runs=[];run=None;step=-1 if r.is_reverse else 1
  for i,v in enumerate(m):
   if v and v['exact']:
    if run and i==run['ros_relative_end0'] and v['uu_pos0']==run['last_uu_pos0']+step:
     run['ros_relative_end0']=i+1;run['last_uu_pos0']=v['uu_pos0']
    else:
     run={'ros_relative_start0':i,'ros_relative_end0':i+1,'first_uu_pos0':v['uu_pos0'],'last_uu_pos0':v['uu_pos0']};runs.append(run)
   else:run=None
  other=[x for x in grouped[key] if not x.is_unmapped and (x.is_secondary or x.is_supplementary)]
  result={'id':key,'ros_chrom':w['chrom'],'ros_start0':int(w['start0']),'uu_chrom':r.reference_name,'uu_start0':r.reference_start,'uu_end0':r.reference_end,'strand':'-' if r.is_reverse else '+','mapq':r.mapping_quality,'cigar':r.cigarstring,'NM':r.get_tag('NM'),'exact_bp':exact,'mismatch_bp':mismatch,'query_unaligned_or_inserted_bp':queryonly,'reference_only_bp':refonly,'additional_reported_alignments':len(other),'exact_runs':runs,'status':'whole_window_exact' if exact==1000 and len(runs)==1 else 'partial_mapping_only_not_window_clearance'}
  assert exact+mismatch+queryonly==1000
  results.append(result);allmaps[key]=m
(O/'UU_base_maps.json').write_text(json.dumps(allmaps));(O/'UU_review.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':results,'limitations':['Primary alignment plus reported secondary/supplementary evidence; not an exhaustive uniqueness proof','Only REF spans fully within one exact run may be considered for simple allele transfer; gaps and mismatches remain explicit','Partial mapping does not satisfy original full-window criterion']},indent=2))
# Diagnose every original chain base, preserving alternative mappings.
by=collections.defaultdict(list);maps={k:[set() for _ in range(1000)] for k in rows}
for k,w in rows.items():by[w['chrom']].append((k,w))
chain=Path('/mnt/d/Jin2024_work/zoonomia_hal/liftover/GCF_014441545.1ToCanFam3.over.chain.gz')
with gzip.open(chain,'rt') as f:
 for line in f:
  x=line.split()
  if not x:continue
  if x[0]=='chain':h=x;tp=int(x[5]);qp=int(x[10]);assert x[4]=='+';continue
  n=int(x[0])
  for key,w in by.get(h[2],[]):
   s=int(w['start0']);e=int(w['end0'])
   for p in range(max(tp,s),min(tp+n,e)):
    q=qp+p-tp;q=q if h[9]=='+' else int(h[8])-1-q;maps[key][p-s].add((h[7],q,h[9]))
  tp+=n;qp+=n
  if len(x)==3:tp+=int(x[1]);qp+=int(x[2])
chainresults=[]
for key,m in maps.items():
 runs=[];last=None
 for i,v in enumerate(m):
  if len(v)!=1:last=None;continue
  chrom,p,strand=next(iter(v));step=1 if strand=='+' else -1
  if last and last['chrom']==chrom and last['strand']==strand and last['ros_relative_end0']==i and last['last_target_pos0']+step==p:last['ros_relative_end0']=i+1;last['last_target_pos0']=p
  else:last={'ros_relative_start0':i,'ros_relative_end0':i+1,'chrom':chrom,'strand':strand,'first_target_pos0':p,'last_target_pos0':p};runs.append(last)
 chainresults.append({'id':key,'unmapped_bp':sum(not v for v in m),'ambiguous_bp':sum(len(v)>1 for v in m),'unique_bp':sum(len(v)==1 for v in m),'unique_contiguous_runs':runs})
(O/'CanFam3_review.json').write_text(json.dumps({'results':chainresults,'limitations':['Chain coordinates only: split runs have not passed reciprocal mapping or sequence identity checks','Unmapped and ambiguous bases remain unresolved']},indent=2))
print(json.dumps([{'id':r['id'],'exact_bp':r['exact_bp'],'mismatch_bp':r['mismatch_bp'],'query_gap_bp':r['query_unaligned_or_inserted_bp'],'reference_gap_bp':r['reference_only_bp'],'other_hits':r['additional_reported_alignments']} for r in results if r['status']!='whole_window_exact'],indent=2))
