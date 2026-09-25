from pathlib import Path
import csv,json,collections
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'06_v1.22.1'
def read(name):return list(csv.DictReader((O/name).open(),delimiter='\t'))
windows=read('review_windows.tsv');mapping={r['window_id']:r for r in read('dog10k_mapping.tsv')};reg={r['window_id']:r for r in read('regulatory_context.tsv')};tiles=read('tile_mappability_summary.tsv');variants=read('small_variant_summary.tsv');records=read('small_variant_records.tsv');raw={};unique=set()
for name in {r['callset'] for r in variants}:
 lines=(O/(name+'.regional.vcf')).read_text().splitlines();raw[name]=[l.split('\t') for l in lines if not l.startswith('#')]
 header=next(l for l in lines if l.startswith('#CHROM'));assert len(header.split('\t'))-9==1987
 for s in [v for v in variants if v['callset']==name]:
  m=mapping[s['window_id']]
  if m['eligible_for_variant_lookup']!='True':assert s['records']=='';continue
  matches=[v for v in raw[name] if v[0]==m['vcf_contig'] and int(v[1])-1<int(m['end']) and int(v[1])-1+len(v[3])>int(m['start'])]
  assert len(matches)==int(s['records']);assert sum(v[6]=='PASS' for v in matches)==int(s['PASS'])
  common=0
  for v in matches:
   info=dict(x.split('=',1) for x in v[7].split(';') if '=' in x);afs=[float(x) for x in info.get('AF','.').split(',') if x!='.'];common+=v[6]=='PASS' and bool(afs) and max(afs)>=.01;unique.add((name,v[0],v[1],v[3],v[4]))
  assert common==int(s['PASS_max_alt_AF_ge_0_01'])
out=[]
for r in windows:
 wid=r['window_id'];m=mapping[wid];r.update({'uu_'+k:v for k,v in m.items() if k not in ['window_id','genes']});r.update({k:v for k,v in reg[wid].items() if k!='window_id'})
 rr=[t for t in tiles if t['window_id']==wid];r['remap_tiles_tested']=sum(int(t['tested']) for t in rr);r['remap_tiles_failed_proxy']=sum(int(t['tested'])-int(t['passed']) for t in rr)
 for v in [v for v in variants if v['window_id']==wid]:
  prefix='SNP' if '.SNPs.' in v['callset'] else 'nonSNP'
  for k in ['records','PASS','PASS_max_alt_AF_ge_0_01']:r[prefix+'_'+k]=v[k]
 r['local_review_decision']='exploratory_only_no_safe_harbor_validation'
 if m['eligible_for_variant_lookup']!='True':r['additional_mapping_flag']='cross_reference_difference_exceeds_existing_NM10_lookup_rule; variants_not_assessed'
 elif r['remap_tiles_failed_proxy']>0:r['additional_mapping_flag']='alternative_tile_alignments_reported; mappability_proxy_not_uniform'
 else:r['additional_mapping_flag']='no_extra_mapping_flag_under_tested_proxy'
 out.append(r)
with (O/'alternatives_integrated_v1221.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
summary=dict(status='passed',windows=len(out),variant_lookup_eligible=sum(m['eligible_for_variant_lookup']=='True' for m in mapping.values()),window_variant_records=len(records),unique_variant_records_in_reviewed_windows=len(unique),tile_tests=sum(r['remap_tiles_tested'] for r in out),tile_failures=sum(r['remap_tiles_failed_proxy'] for r in out),windows_with_tile_failures=[r['window_id'] for r in out if r['remap_tiles_failed_proxy']],all_regulatory_overlap_zero=all(all(int(r[k])==0 for k in ['epic_promoter_enhancer_bp','lncrna_smallrna_overlap_bp','mirna_overlap_bp','external_regulatory_overlap_bp']) for r in out),validation=['Raw VCF text parser reproduces overlapping REF-span record counts, PASS counts and AF>=1% counts independently of pysam object extraction','1987 genotype columns in both source VCFs','All window records REF checked against UU reference during extraction','Variant summaries keyed by unique window ID; ineligible lookup remains missing','Mapping, regulatory and variant tables integrated without promoting any safe-harbor status'],limitations=['Independently parsing cached records verifies arithmetic, not remote index completeness','Reference-span overlap convention, including anchors and reference-consuming indels','Tile tests include duplicates in overlapping windows; not independent loci'])
(O/'integration_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(json.dumps([{k:r[k] for k in ['window_id','genes','start','end','SNP_PASS','nonSNP_PASS','SNP_PASS_max_alt_AF_ge_0_01','nonSNP_PASS_max_alt_AF_ge_0_01','remap_tiles_failed_proxy']} for r in out],indent=2))
