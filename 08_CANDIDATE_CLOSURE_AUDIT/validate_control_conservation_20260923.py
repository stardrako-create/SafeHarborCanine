from pathlib import Path
import json,hashlib,collections,datetime
A=Path(__file__).resolve().parent;O=A/'CONTROL_LOCUS_REVIEW_20260922/conservation_20260923'
s=json.loads((O/'status.json').read_text());results=[]
for path in sorted(O.glob('*.maf')):
 key=path.stem;species=set();support={};block=[]
 def emit():
  seq=[l.split() for l in block if l.startswith('s\t') or l.startswith('s ')];dogs=[r for r in seq if r[1].startswith('Canis_lupus_familiaris.')]
  if not dogs:return
  assert len(dogs)==1;dog=dogs[0];assert dog[4]=='+';p=int(dog[2])
  for r in seq:species.add(r[1].split('.')[0]);assert len(r[6])==len(dog[6])
  for i,b in enumerate(dog[6]):
   if b=='-':continue
   assert p not in support
   support[p]=len({r[1].split('.')[0] for r in seq if not r[1].startswith('Canis_lupus_familiaris.') and r[6][i].upper() in 'ACGT'})
   p+=1
 for line in path.read_text().splitlines()+['']:
  if not line.strip():emit();block=[]
  else:block.append(line)
 assert len(support)==1000
 completed=next((r for r in s['completed'] if r['id']==key),None)
 wig=O/(key+'.wig')
 if completed:
  assert hashlib.sha256(path.read_bytes()).hexdigest()==completed['maf_sha256'];assert hashlib.sha256(wig.read_bytes()).hexdigest()==completed['wig_sha256']
 result={'id':key,'reference_bp':len(support),'species_total':len(species),'bp_with_non_dog_ACGT':sum(v>0 for v in support.values()),'bp_with_at_least_10_non_dog_species':sum(v>=10 for v in support.values()),'status':'scores_available_with_alignment_support_reported' if completed else 'unresolved'}
 if not completed and species=={'Canis_lupus_familiaris'} and not wig.read_text().strip():result['status']='no_comparative_alignment_no_conservation_score'
 results.append(result)
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':results,'note':'Original execution status and partial outputs preserved. bg1k_0117 has only reference sequence; empty phyloP output is missing comparative evidence, not zero conservation. Alignment support counts are descriptive, not new pass thresholds.'}
(O/'alignment_support_review.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(results,indent=2))
