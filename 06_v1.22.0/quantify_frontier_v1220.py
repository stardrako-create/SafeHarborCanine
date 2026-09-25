from pathlib import Path
import csv,json,struct,pysam
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.0';OLD=P/'06_v1.21.7'
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
assign=read(P/'06_v1.21.9/RNA_cluster_assignments.tsv');cols=[c for c in assign[0] if c.endswith('_T_enriched')];groups=[{r['barcode'] for r in assign if r[c]=='1'} for c in cols]
masks={'T_consensus':set.intersection(*groups),'T_union':set.union(*groups),'all_QC':{r['barcode'] for r in assign}}
queries=json.loads((O/'regional_fragment_queries.json').read_text());fragments={};valid=[];magic=b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00'
for q in queries['queries']:
 d=(OLD/'range_cache'/f"{q['range_start']}-{q['range_end']}.bgzfpart").read_bytes();a=d.find(magic);assert a>=0;b=a
 while b+18<=len(d) and d[b:b+16]==magic:
  n=struct.unpack_from('<H',d,b+16)[0]+1
  if b+n>len(d):break
  b+=n
 p=O/(q['target']['name']+'.verified_chunk.bgz')
 with pysam.BGZFile(str(p),'wb'):pass
 eof=p.read_bytes();p.write_bytes(d[a:b]+eof)
 with pysam.BGZFile(str(p),'rb') as f:lines=f.read().decode().split('\n')[1:-1]
 rows=[l.split('\t') for l in lines if l and not l.startswith('#')];t=q['target'];hits=[r for r in rows if r[0]==t['chrom'] and int(r[1])<t['end'] and int(r[2])>t['start']]
 saved=[l.split('\t') for l in (O/q['raw_file']).read_text().splitlines() if l];assert hits==saved;assert len({tuple(r[:4]) for r in saved})==len(saved)
 fragments[t['name'].removesuffix('_frontier_span')]=saved;valid.append(dict(target=t['name'],HTSlib_match=True,records=len(saved)))
rows=read(O/'frontier_canfam6_mapping.tsv');summ=[]
blockmaps={(r['genes'],r['start'],r['end']):r for r in json.loads((O/'frontier_shared_mapping_blocks.json').read_text())};unique_hits={}
for r in rows:
 for m in masks:r[m+'_fragment_starts']='';r[m+'_nuclei_with_start']=''
 r['T_review_status']='not_quantified_mapping_unresolved'
 bm=blockmaps[(r['genes'],int(r['start']),int(r['end']))];r['shared_mapped_bp']=bm['shared_unique_destination_bp'];r['source_mapping_agreement_fraction']=bm['source_bp_with_two_route_agreement']/bm['source_bp'];r['shared_blocks_count']=len(bm['blocks'])
 data=fragments[r['genes'].split('/')[0]];coords={p for c,a,b in bm['blocks'] for p in range(a,b)};assert len(coords)==bm['shared_unique_destination_bp']
 for m,barcodes in masks.items():
  h=[v for v in data if v[3] in barcodes and int(v[1]) in coords];r[m+'_shared_fragment_starts']=len(h);r[m+'_shared_nuclei_with_start']=len({v[3] for v in h})
  if m=='T_consensus':unique_hits.setdefault(r['genes'],set()).update(tuple(v[:4]) for v in h)
 if r['canfam6_mapping_status']=='two_routes_unique_identical_contiguous':
  data=fragments[r['genes'].split('/')[0]];s=int(r['canfam6_start']);e=int(r['canfam6_end'])
  for m,barcodes in masks.items():
   h=[v for v in data if v[3] in barcodes and s<=int(v[1])<e];r[m+'_fragment_starts']=len(h);r[m+'_nuclei_with_start']=len({v[3] for v in h})
  r['T_review_status']='exploratory_sparse_signal_requires_replication' if r['T_consensus_fragment_starts'] else 'no_start_observed_not_proof_closed'
for g in sorted({r['genes'] for r in rows}):
 rr=[r for r in rows if r['genes']==g];ok=[r for r in rr if r['T_consensus_fragment_starts']!=''];positive=[r for r in ok if r['T_consensus_fragment_starts']>0]
 summ.append(dict(genes=g,frontier_windows=len(rr),quantified_full_contiguous=len(ok),partial_or_discontinuous=len(rr)-len(ok),T_positive_full_windows=len(positive),T_positive_shared_windows=sum(r['T_consensus_shared_fragment_starts']>0 for r in rr),min_source_mapping_agreement=min(r['source_mapping_agreement_fraction'] for r in rr),max_T_shared_starts=max(r['T_consensus_shared_fragment_starts'] for r in rr),max_T_shared_nuclei=max(r['T_consensus_shared_nuclei_with_start'] for r in rr),unique_T_fragment_records_across_windows=len(unique_hits[g]),unique_T_nuclei_across_windows=len({v[3] for v in unique_hits[g]})))
def write(p,r):
 with p.open('w') as f:w=csv.DictWriter(f,fieldnames=list(r[0]),delimiter='\t');w.writeheader();w.writerows(r)
write(O/'frontier_T_evidence_v1220.tsv',rows)
result=dict(validation=valid,summary=summ,masks={k:len(v) for k,v in masks.items()},limitations=['37 of 64 windows meet exact contiguous two-route mapping; shared-base counts also provided for all64, with coverage reported; missing bases not assigned zero','Shared-base count is a different interval definition from full1kb count','Start-based counts; not all overlapping fragments','Overlapping windows and nuclei from one tumor donor are not independent biological replicates','No full-depth normalization or significance claim','T groups exploratory RNA-based identity; source ATAC QC still incomplete'])
(O/'frontier_T_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
