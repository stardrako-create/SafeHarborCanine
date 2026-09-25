from pathlib import Path
import csv,json,collections,struct,pysam,io
O=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino/06_v1.21.7')
cells={r['barcode']:r for r in csv.DictReader((O/'cell_marker_masks.tsv').open(),delimiter='\t')};masks=['all_QC','T_marker_2plus','T_marker_strict','T_marker_broad'];result=[];validation=[]
q=json.loads((O/'regional_fragment_queries.json').read_text());MAGIC=b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00'
for item in q['queries']:
 t=item['target'];d=(O/'range_cache'/f"{item['range_start']}-{item['range_end']}.bgzfpart").read_bytes();a=d.find(MAGIC);p=a
 while p+18<=len(d) and d[p:p+16]==MAGIC:
  size=struct.unpack_from('<H',d,p+16)[0]+1
  if p+size>len(d):break
  p+=size
 archive=O/(t['name'].replace('/','_')+'.verified_chunk.bgz')
 with pysam.BGZFile(str(archive),'wb') as out:
  # Recompression is only a container step; HTSlib independently decodes the original complete blocks below.
  pass
 eof=archive.read_bytes();archive.write_bytes(d[a:p]+eof)
 with pysam.BGZFile(str(archive),'rb') as inp:lines=inp.read().decode().split('\n')[1:-1]
 rows=[l.split('\t') for l in lines if l and not l.startswith('#')];selected=[r for r in rows if r[0]==t['chrom'] and int(r[1])<t['end'] and int(r[2])>t['start']]
 saved=[l.split('\t') for l in (O/item['raw_file']).read_text().splitlines() if l];assert selected==saved
 validation.append(dict(target=t['name'],HTSlib_decoded_rows_match=True,queried_overlap_records=len(selected)))
 for mask in masks:
  eligible={b for b,c in cells.items() if c[mask]=='1'};r=[v for v in selected if v[3] in eligible];starts=[v for v in r if t['start']<=int(v[1])<t['end']];barcodes={v[3] for v in starts}
  result.append(dict(target=t['name'],kind=t['kind'],mask=mask,nuclei=len(eligible),unique_fragment_records_starting_in_window=len(starts),nuclei_with_fragment_start=len(barcodes),unique_fragment_records_overlapping_fetched_window=len(r),read_support_sum_starting_in_window=sum(int(v[4]) for v in starts),fraction_nuclei_with_start=len(barcodes)/len(eligible),ATAC_matrix_peak_count_sum=sum(int(cells[b]['ATAC_peak_counts']) for b in eligible)))
with (O/'local_T_marker_fragment_counts.tsv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');wr.writeheader();wr.writerows(result)
(O/'fragment_quantification_validation.json').write_text(json.dumps(dict(queries=validation,unique_records_not_PCR_support_used_as_primary=True,limitations=['Coordinate-sorted BGZF searched by byte probes; only regional data downloaded','Complete start-based counts within bracketed intervals; overlap counts could omit extremely long fragments from before fetched range','RNA masks are exploratory, not validated T-cell clusters','Marker controls are feature-associated intervals, not confirmed promoters','One tumor donor; nuclei do not supply biological replication','No normalization to full genomic fragment depth; matrix peak count sum is a different quantity']),indent=2))
print(json.dumps([r for r in result if r['mask']=='T_marker_strict'],indent=2))
