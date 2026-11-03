from pathlib import Path
import json,pysam,hashlib,collections,datetime
O=Path(__file__).parent
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
D=O/('allele_projection_'+stamp);D.mkdir(exist_ok=False)
blocks={r['id']:r for r in json.loads((O/'basewise_projection.json').read_text())['results']}
rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
results=[];sources=[]
with pysam.FastaFile(str(O/'UU_projection/controls_ROS.fa')) as ros:
 for folder in ['variant_catalogue_query_drop_samples','variant_catalogue_query_remaining3']:
  for p in sorted((O/folder).glob('*__*.json')):
   data=p.read_bytes();j=json.loads(data);w=blocks[j['window']];seq=ros.fetch(j['window']).upper()
   sources.append(dict(path=str(p),sha256=hashlib.sha256(data).hexdigest()))
   for a in j['alleles']:
    s=a['pos1']-1;e=s+len(a['ref']);candidate=[b for b in w['exact_blocks'] if b['UU_start0']<=s<e<=b['UU_end0']]
    result=dict(window=j['window'],source=j['source'],allele=a,status='needs_boundary_or_reference_reconciliation')
    if a['reference_verified'] and a['sequence_allele'] and len(candidate)==1:
     b=candidate[0];off=b['ROS_offset_start0']+s-b['UU_start0'] if w['strand']=='+' else b['ROS_offset_end0']-(e-b['UU_start0'])
     ref=a['ref'].upper() if w['strand']=='+' else rc(a['ref'].upper());alt=a['alt'].upper() if w['strand']=='+' else rc(a['alt'].upper())
     assert seq[off:off+len(ref)]==ref
     result.update(status='exact_block_REF_verified_normalization_pending',ROS_offset0=off,ROS_ref=ref,ROS_alt=alt)
    results.append(result)
(D/'review.json').write_text(json.dumps(dict(utc=stamp,sources=sources,alleles=results,counts=dict(collections.Counter(x['status'] for x in results)),limitations=['Snapshot of available query receipts; missing queries not negative','Projected alleles not normalized yet','No individual genotype or guide specificity conclusion']),indent=2))
print(json.dumps({'directory':str(D),'queries':len(sources),'alleles':len(results),'counts':dict(collections.Counter(x['status'] for x in results))}))
