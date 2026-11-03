from pathlib import Path
import json
A=Path(__file__).resolve().parent;O=A/'W11_LOCAL_MODELS_20261002';hits=json.loads((O/'candidate_pair_hits.json').read_text());rows=[]
for name in sorted({k.rsplit('_mate',1)[0] for k in hits}):
 for model in json.loads((O/'models.json').read_text()):
  hh=[]
  for m in [1,2]:
   group=[x for x in hits.get(name+'_mate'+str(m),[]) if x['model']==model]
   if not group:break
   hh.append(max(group,key=lambda x:x['AS']))
  if len(hh)==2:
   left,right=sorted(hh,key=lambda x:x['start0']);rows.append({'name':name,'model':model,'sum_AS':sum(x['AS'] for x in hh),'inward_FR':not left['reverse'] and right['reverse'],'outer_span':max(x['end0'] for x in hh)-min(x['start0'] for x in hh),'mates':hh})
(O/'pair_geometry.json').write_text(json.dumps({'rows':rows,'limits':['Local model geometry only, no calibrated insert size or genomic competition','Best reported alignment per mate/model; ties may remain','Does not demonstrate second junction']},indent=2))
for name in sorted({x['name'] for x in rows}):
 best=max(x['sum_AS'] for x in rows if x['name']==name);rr=[x for x in rows if x['name']==name and x['sum_AS']==best];print(name,best,'best_models',len(rr),'inward',sum(x['inward_FR'] for x in rr))