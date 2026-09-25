from pathlib import Path
import csv,json,collections,bisect,random,math,pysam,hashlib
P=Path(__file__).parent.parent;O=Path(__file__).parent/'CONTROL_LOCUS_REVIEW_20260922';O.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def bed(p):
 d=collections.defaultdict(list)
 for l in p.open():
  if l.startswith('#') or not l.strip():continue
  x=l.rstrip().split('\t');d[x[0]].append((int(x[1]),int(x[2]),x[3:]))
 return d
def dist(items,s,e):return min((max(a-e,s-b,0) for a,b,*_ in items),default=None)
def overlap(items,s,e):return any(a<e and b>s for a,b,*_ in items)
genes=bed(P/'05_SHIP/canine_all_genes.bed');coding=bed(P/'05_SHIP/canine_all_genes_stranded.bed');mir=bed(P/'05_SHIP/canine_miRNA.bed');lnc=bed(P/'05_SHIP/canine_lncRNA_smallRNA.bed');reg=bed(P/'05_SHIP/ehsan_regulatory_elements_ROS.bed');peaks=bed(P/'06_v1.21.5/consensus_peaks_76_min39_gap75.bed')
risk={r['gene_symbol'] for r in read(P/'05_SHIP/canine_risk_genes.tsv')};riskbed={c:[r for r in v if r[2] and r[2][0] in risk] for c,v in genes.items()}
tss={c:[(a if rest[-1]=='+' else b-1,(a if rest[-1]=='+' else b-1)+1,[]) for a,b,rest in v] for c,v in coding.items()}
bg=[r for r in read(P/'06_v1.21.8/local_background_windows.tsv') if r['width']=='1000'];vals=sorted(float(r['value']) for r in bg if r['value'] not in ('','None') and math.isfinite(float(r['value'])))
ref=pysam.FastaFile(str(P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'))
wseq=ref.fetch('NC_051811.1',17255706,17256706).upper();wgc=(wseq.count('G')+wseq.count('C'))/len(wseq)
rows=[]
for i,r in enumerate(bg):
 c=r['chrom'];s=int(r['start']);e=int(r['end']);assert e-s==1000
 seq=ref.fetch(c,s,e).upper();canonical=len(seq)==1000 and not(set(seq)-set('ACGT'));gc=(seq.count('G')+seq.count('C'))/len(seq)
 observed=r['value'] not in ('','None') and math.isfinite(float(r['value']));pct=100*bisect.bisect_right(vals,float(r['value']))/len(vals) if observed else None
 td=dist(tss.get(c,[]),s,e);md=dist(mir.get(c,[]),s,e);rd=dist(riskbed.get(c,[]),s,e)
 flags=[]
 for key,yes in [('missing_ATAC',not observed),('noncanonical_reference',not canonical),('gene_overlap',overlap(genes.get(c,[]),s,e)),('coding_5prime_under_50kb',td is not None and td<50000),('miRNA_under_300kb',md is not None and md<300000),('risk_gene_under_300kb',rd is not None and rd<300000),('lnc_smallRNA_overlap',overlap(lnc.get(c,[]),s,e)),('regulatory_overlap',overlap(reg.get(c,[]),s,e)),('consensus_peak_overlap',overlap(peaks.get(c,[]),s,e)),('annotation_contig_missing',c not in genes or c not in coding)]:
  if yes:flags.append(key)
 rows.append(dict(id=f'bg1k_{i:04d}',chrom=c,start0=s,end0=e,atac_mean=r['value'],atac_percentile=pct,GC=gc,GC_distance_w01=abs(gc-wgc),coding_TSS_distance=td,miRNA_distance=md,risk_gene_distance=rd,screen_flags=';'.join(flags),stage='preliminary_annotation_screen_only',pending='alternative_transcripts;TAD;EpiC;repeats;conservation;variants;sequence_specificity;donor_genotype'))
valid=[r for r in rows if not r['screen_flags']];rng=random.Random(20260922);randomset=rng.sample(valid,min(10,len(valid)))
low=sorted([r for r in valid if r['atac_percentile']<=20],key=lambda r:(r['GC_distance_w01'],r['id']))[:10]
def write(name,rs,fields):
 with (O/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rs)
write('all_3000_windows_annotation_screen.tsv',rows,list(rows[0]));write('random_preliminary_review_queue.tsv',randomset,list(rows[0]));write('low_ATAC_preliminary_review_queue.tsv',low,list(rows[0]))
summary=dict(background_windows=len(rows),observed_ATAC=len(vals),without_flags_in_limited_screen=len(valid),random_queue=len(randomset),low_ATAC_queue=len(low),shared_ids=sorted({r['id'] for r in randomset}&{r['id'] for r in low}),w01_GC=wgc,random_seed=20260922,low_ATAC_rule='<=20th percentile of observed 1kb background; exploratory comparator only, not candidate eligibility threshold',decision='Review queues only; no locus released for editing',limitations=['Random sample conditional on limited annotation filters; not unrestricted genomic random sample','Low-ATAC ranked by GC distance, not fully covariate-matched','Gene-level 5prime screen does not complete alternative-transcript review','Remaining mandatory evidence explicitly pending; no missing layer treated as passing'])
(O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
