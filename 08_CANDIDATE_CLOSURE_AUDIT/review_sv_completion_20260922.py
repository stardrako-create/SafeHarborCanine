from pathlib import Path
import json, hashlib, pysam, datetime
A=Path(__file__).parent/'SV_PBGV000010/alignment_chunked'
receipts=[]
for p in sorted(A.glob('chunk_*/validated.json')):
 r=json.loads(p.read_text()); h=hashlib.sha256()
 with (p.parent/'retained.bam').open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
 assert h.hexdigest()==r['bam_sha256'],str(p)
 assert r['primary']==2*r['pairs'] and r['read1']==r['read2']==r['pairs']
 receipts.append(r)
assert len(receipts)==47 and sum(r['pairs'] for r in receipts)==185353881
assert [r['start_pair0'] for r in receipts]==list(range(0,185353881,4000000))
bam=A/'selected_markdup.bam'
pysam.quickcheck(str(bam))
chrom='NC_049253.1'; boundaries=[4707938,33613308]; results=[]
with pysam.AlignmentFile(str(bam),'rb') as f:
 assert f.has_index()
 for i,pos in enumerate(boundaries):
  other=boundaries[1-i]; paired=set(); split=set(); eligible=set()
  for r in f.fetch(chrom,pos-10000,pos+10000):
   if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate or r.mapping_quality<30:continue
   eligible.add(r.query_name)
   if not r.mate_is_unmapped and r.next_reference_name==chrom and abs(r.next_reference_start-other)<=10000:paired.add(r.query_name)
   if r.has_tag('SA'):
    for sa in r.get_tag('SA').rstrip(';').split(';'):
     c,p,s,cigar,mq,nm=sa.split(',')
     if c==chrom and abs(int(p)-1-other)<=10000 and int(mq)>=30:split.add(r.query_name)
  results.append(dict(boundary0=pos,eligible_unique_names=len(eligible),mate_links_to_other_boundary=sorted(paired),SA_links_to_other_boundary=sorted(split)))
result=dict(time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),validated_receipts=47,pairs=185353881,all_chunk_hashes_match=True,final_bam_quickcheck=True,final_bam_index=True,boundary_link_screen=results,limitations=['Screen uses +/-10 kb around catalogue coordinates and MAPQ30 nonduplicate primary anchors','SA links use supplementary alignment start, not reconstructed junction coordinates','Mate mapping quality and orientation not used to establish SV support','Absence of links in this screen does not exclude shifted breakpoints or complex/mosaic events','Coverage GC/mappability matching and control review remain pending'])
(A/'receipt_and_boundary_link_review_20260922.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
