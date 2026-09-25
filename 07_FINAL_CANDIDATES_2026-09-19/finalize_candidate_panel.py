from pathlib import Path
import csv,json,hashlib,shutil,datetime
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'07_FINAL_CANDIDATES_2026-09-19';O.mkdir(exist_ok=True)
def read(p):return list(csv.DictReader((P/p).open(),delimiter='\t'))
rows=read('06_v1.22.2/alternatives_integrated_v1222.tsv');assert len(rows)==16 and len({r['window_id'] for r in rows})==16
sources=['06_v1.22.2/alternatives_integrated_v1222.tsv','06_v1.22.1/small_variant_summary.tsv','06_v1.22.1/tile_mappability_summary.tsv','06_v1.22.1/regulatory_context.tsv','06_v1.22.1/dog10k_mapping.tsv','06_v1.22.0/frontier_T_evidence_v1220.tsv','06_v1.22.3/depth_matched_local_counts.tsv','06_v1.22.9/all40_overview.tsv','06_v1.21.8/CURRENT.md','06_v1.22.2/LOC_mapping_resolution.json','06_v1.22.2/LOC_supplementary_lookup.json','06_v1.21.9/NPNT_span_review.json']
reg={r['window_id']:r for r in read(sources[3])};maps={r['window_id']:r for r in read(sources[4])};tiles=read(sources[2]);variants=read(sources[1]);front=read(sources[5]);depth=read(sources[6])
assert len(front)==64
for r in rows:
 k=r['window_id'];assert int(r['end'])-int(r['start'])==1000
 for f in ['nearest_gene_body','gene_body_gap_bp','nearest_gene_5prime','gene_5prime_gap_bp','epic_promoter_enhancer_bp','lncrna_smallrna_overlap_bp','mirna_overlap_bp','external_regulatory_overlap_bp']:assert r[f]==reg[k][f],(k,f)
 t=[x for x in tiles if x['window_id']==k];assert sum(int(x['tested']) for x in t)==int(r['remap_tiles_tested']);assert sum(int(x['tested'])-int(x['passed']) for x in t)==int(r['remap_tiles_failed_proxy'])
 if k!='w00':
  for c,prefix in [('AutoAndXPAR.SNPs.vqsr99.vcf.gz','SNP'),('AutoAndXPAR.nonSNPs.filter.vcf.gz','nonSNP')]:
   v=next(x for x in variants if x['window_id']==k and x['callset']==c)
   for f in ['PASS','records','PASS_max_alt_AF_ge_0_01']:assert v[f]==r[prefix+'_'+f],(k,f)
  assert r['uu_sv_records']==maps[k]['sv_records']
 f=next(x for x in front if x['chrom']==r['chrom'] and x['start']==r['start'] and x['end']==r['end'])
 for col in ['T_consensus_shared_fragment_starts','T_union_shared_fragment_starts','shared_mapped_bp']:assert f[col]==r[col]
decisions={
'w00':('NOT_ADVANCED','Only 5% RRBS coverage, max one dog; 204 bp repeat; cross-assembly core differs and flanks split; 923/1000 bp shared for T evidence.'),
'w01':('PRIORITY_1_VALIDATION','No tested repeat overlap; full cross-reference mapping; no catalogue SV hit; better RRBS coverage than w11; two sparse T starts. Methylation remains high and CAR-T accessibility unproven.'),
'w02':('SAME_LOCUS_BACKUP','Overlaps w01 by 500 bp; not an independent locus. Higher bulk ATAC but more variant records, 400 bp nonquiescent annotation and only one T start.'),
'w03':('RESERVE_NOT_FIRST_LINE','One T start on 998 shared bp; 208 bp repeats; RRBS methylation 88.27% on 10% coverage; no catalogue SV hit. Not superior overall despite bulk ATAC percentile 90.95.'),
'w04':('NOT_ADVANCED','Zero consensus T starts; one union-only start. RRBS max three dogs, repeats and nonquiescent annotation.'),
'w05':('NOT_ADVANCED','Zero consensus T starts; one union-only start. RRBS max three dogs; near NTF3 body; nonquiescent annotation.'),
'w06':('NOT_ADVANCED','Zero consensus T starts; one union-only start. Zero methylation measured on only 10% window with max one dog; not robust hypomethylation.'),
'w07':('NOT_ADVANCED','Two SV catalogue hits including smaller deletion; 996 shared bp; high observed methylation and seven PASS nonSNP records.'),
'w08':('NOT_ADVANCED','Same NPNT SV context; high observed methylation (92.68%) and repeat overlap.'),
'w09':('NOT_ADVANCED','Overlapping NPNT alternative; 494 bp repeats and 81.96% observed methylation; not preferred over w11.'),
'w10':('SAME_LOCUS_BACKUP','Overlaps w11 by 750 bp; same SV context; 363 bp repeats vs157 and one T start vs2; tiny bulk ATAC advantage not decisive.'),
'w11':('PRIORITY_2_CONDITIONAL','Full reference mapping; 103/103 tile checks; two sparse T starts; bulk ATAC percentile89.50. Conditional on resolving NPNT SV context; RRBS coverage only10%, observed methylation63.93%.'),
'w12':('SAME_LOCUS_BACKUP','Overlaps w11 by750 bp; same SV context; more repeats/nonSNPs and one T start; small ATAC difference not decisive.'),
'w13':('NOT_ADVANCED','Eight failed conservative mapping tiles, 424 bp repeats,96.20% observed methylation; shared NPNT SV concern.'),
'w14':('NOT_ADVANCED','Four T starts among490 all-QC starts do not demonstrate T enrichment; RRBS coverage3.9%, max2 dogs,100% observed methylation; SV concern.'),
'w15':('NOT_ADVANCED','Three T starts among342 all-QC starts; RRBS5%, max1 dog,100% observed methylation; SV concern.')}
def write(path,records):
 with path.open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
out=[]
for r in rows:out.append(dict(window_id=r['window_id'],decision=decisions[r['window_id']][0],reason=decisions[r['window_id']][1],**{k:v for k,v in r.items() if k!='window_id'}))
write(O/'all16_decisions_and_evidence.tsv',out)
selected=[next(r for r in out if r['window_id']==k) for k in ['w01','w11']];write(O/'shortlist_evidence.tsv',selected)
reserve=next(r for r in out if r['window_id']=='w03');write(O/'reserve_evidence.tsv',[reserve])
for filename,data in [('shortlist_ROS_0based.bed',selected),('reserve_ROS_0based.bed',[reserve])]:
 (O/filename).write_text(''.join(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['window_id']}\n" for r in data))
cross=[]
for r in selected+[reserve]:
 cross.append(dict(window_id=r['window_id'],ROS_accession=r['chrom'],ROS_start0=int(r['start']),ROS_end0=int(r['end']),ROS_start1=int(r['start'])+1,ROS_end1=int(r['end']),UU_contig=r['uu_uu_contig'],UU_start0=r['uu_start'],UU_end0=r['uu_end'],UU_strand_relative_ROS=r['uu_strand'],canFam6_contig=r['canfam6_chrom'],canFam6_start0=r['canfam6_start'],canFam6_end0=r['canfam6_end'],shared_T_mapping_bp=r['shared_mapped_bp']))
write(O/'coordinate_crosswalk.tsv',cross)
remaining=[r for r in front if not any(r['chrom']==x['chrom'] and r['start']==x['start'] and r['end']==x['end'] for x in rows)];assert len(remaining)==48
write(O/'remaining48_no_new_biological_veto.tsv',[dict(review_disposition='Not promoted in current evidence-bounded prioritization; low/zero T signal is NOT a veto',**r) for r in remaining])
(O/'input_manifest_sha256.json').write_text(json.dumps({s:hashlib.sha256((P/s).read_bytes()).hexdigest() for s in sources},indent=2))
(O/'integration_checks.json').write_text(json.dumps(dict(status='passed',alternatives=16,prior_frontier=64,unpromoted_other_windows=48,selected=['w01','w11'],reserve=['w03'],scope='Decision panel for experimental validation, not proven safe harbors or guide-ready insertion sites',source_fields_crosschecked=True),indent=2))
print('Source checks passed; provisional panel assembled. Sequence and report checks pending.')
