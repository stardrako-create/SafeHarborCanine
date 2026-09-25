from pathlib import Path
import csv,json,hashlib,shutil
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'08_CANDIDATE_CLOSURE_AUDIT'
read=lambda p:list(csv.DictReader(p.open(encoding='utf-8'),delimiter='\t'))
segments=[r for r in read(O/'all_short_segments_exact_uniqueness.tsv') if int(r['length'])==20]
variants=[r for r in read(O/'normalized_variant_alleles_with_AF_ROS.tsv') if r['filter']=='PASS']
prior=read(O/'twenty_base_sequence_constraints.tsv');assert len(segments)==len(prior)==1962
bykey={(r['window_id'],r['start0'],r['end0']):r for r in prior}
summ=[]
for wid in ['w01','w11']:
 counts=dict(window=wid,tested=0,nonunique=0,variant_overlap=0,unflagged=0)
 for r in segments:
  if r['window_id']!=wid:continue
  hits=sorted({v['id'] for v in variants if v['window_id']==wid and v['chrom']==r['chrom'] and int(v['start0'])<int(r['end0']) and int(v['end0'])>int(r['start0'])})
  n=int(r['exact_locus_count_both_strands']);old=bykey[(wid,r['start0'],r['end0'])]
  assert int(old['exact_locus_count'])==n and int(old['PASS_variant_records'])==len(hits)
  assert set(old['PASS_variant_record_ids'].split(';'))-set([''])==set(hits)
  good=n==1 and not hits
  assert (old['unflagged_by_two_limited_checks']=='True')==good
  counts['tested']+=1;counts['nonunique']+=n>1;counts['variant_overlap']+=bool(hits);counts['unflagged']+=good
 summ.append(counts)
sv=json.loads((P/'06_v1.21.6/NPNT_SV_carriers_review.json').read_text())
ids={r['sample'] for r in sv['carriers']}
meta=read(P/'06_v1.21.6/dog10K-alignment-sample-table.2022-02-23-v8.txt.UPDATESRA.txt')
carriers=[{k:r[k] for k in ['sampleName','Breed/Type','BIOSAMID','runID','effectiveAutosomalMeanCoverage']} for r in meta if r['sampleName'] in ids]
assert len(carriers)==4
alignments=[str(p.relative_to(P)) for p in P.rglob('*') if p.suffix.lower() in ['.bam','.cram']]
result=dict(status='completed',all_1962_rows_recomputed_identically=True,summary=summ,SV_carrier_raw_data_accessions_from_local_metadata=carriers,local_project_bam_cram=alignments,scope='Project tree only; public run identifiers are metadata, not verified availability or downloaded evidence',limitations=['No nuclease selection found in project markdown/yaml/json keyword search','No prospective donor genotype identified; catalogue carriers are not the intended CAR-T donors','No raw SV read analysis performed'])
(O/'constraints_reproduction_and_input_inventory.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(__file__,O/Path(__file__).name)
print(json.dumps(result,indent=2))
