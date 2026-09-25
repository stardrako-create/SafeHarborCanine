from pathlib import Path
import urllib.request,gzip,struct,csv,json,hashlib
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino');O=P/'06_v1.21.6'
base='https://kiddlabshare.med.umich.edu/dog10K/SNP_and_indel_calls_2021-10-17/'
windows=list(csv.DictReader((P/'06_v1.21.5/dog10k_local_mapping.tsv').open(),delimiter='\t'))
def get(url,a=None,b=None):
 req=urllib.request.Request(url,headers={} if a is None else {'Range':f'bytes={a}-{b}'})
 with urllib.request.urlopen(req,timeout=45) as r:
  if a is not None:
   assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {a}-{b}/'),dict(r.headers)
  data=r.read(25000001);assert len(data)<=25000000
  return data
def bgzf(data):
 pos=0;out=[]
 while pos+18<=len(data):
  assert data[pos:pos+2]==b'\x1f\x8b'
  xl=struct.unpack_from('<H',data,pos+10)[0];extra=data[pos+12:pos+12+xl];i=0;size=None
  while i<len(extra):
   tag=extra[i:i+2];n=struct.unpack_from('<H',extra,i+2)[0]
   if tag==b'BC':size=struct.unpack_from('<H',extra,i+4)[0]+1
   i+=4+n
  assert size
  if pos+size>len(data):break
  out.append((pos,gzip.decompress(data[pos:pos+size])));pos+=size
 return out
def readindex(data):
 d=gzip.decompress(data);assert d[:4]==b'TBI\1';off=4
 def unpack(fmt):
  nonlocal off
  out=struct.unpack_from('<'+fmt,d,off);off+=struct.calcsize('<'+fmt);return out
 nr,fmt,cseq,cbeg,cend,meta,skip,ln=unpack('8i');names=d[off:off+ln].rstrip(b'\0').decode().split('\0');off+=ln;refs={}
 assert fmt==2
 for name in names:
  bins={}
  for _ in range(unpack('i')[0]):
   bid,nc=unpack('Ii');bins[bid]=[unpack('QQ') for _ in range(nc)]
  ni=unpack('i')[0];linear=list(unpack('Q'*ni)) if ni else []
  refs[name]=(bins,linear)
 return refs
def binsfor(a,b):
 b-=1;ans=[0]
 for off,shift in [(1,26),(9,23),(73,20),(585,17),(4681,14)]:ans.extend(range(off+(a>>shift),off+(b>>shift)+1))
 return ans
allsummary=[]
for name in ['AutoAndXPAR.SNPs.vqsr99.vcf.gz','AutoAndXPAR.nonSNPs.filter.vcf.gz']:
 print('Index',name,flush=True);idx=get(base+name+'.tbi');(O/(name+'.tbi')).write_bytes(idx);refs=readindex(idx)
 headerdata=get(base+name,0,1048575);header=b''.join(x[1] for x in bgzf(headerdata)).decode();headers=[l for l in header.splitlines() if l.startswith('#')];assert any(l.startswith('#CHROM') for l in headers)
 (O/(name+'.header.txt')).write_text('\n'.join(headers)+'\n');nsamples=len(next(l for l in headers if l.startswith('#CHROM')).split('\t'))-9
 out=[];raw=[];requests=[]
 for w in windows:
  chrom=w['vcf_contig'];a=int(w['start']);b=int(w['end']);bins,linear=refs[chrom];minoff=linear[min(a>>14,len(linear)-1)] if linear else 0
  chunks=sorted(set((s,e) for bid in binsfor(a,b) for s,e in bins.get(bid,[]) if e>minoff));merged=[]
  for s,e in chunks:
   if merged and s<=merged[-1][1]:merged[-1][1]=max(e,merged[-1][1])
   else:merged.append([s,e])
  found={}
  for s,e in merged:
   lo=s>>16;hi=(e>>16)+65535;assert hi-lo<25000000
   data=get(base+name,lo,hi);requests.append(dict(chrom=chrom,start_byte=lo,end_byte=hi,sha256=hashlib.sha256(data).hexdigest()))
   blocks=bgzf(data);text=b''.join(block for _,block in blocks)[s&65535:].decode()
   for l in text.split('\n')[:-1]:
    v=l.split('\t')
    if len(v)<8 or v[0]!=chrom:continue
    start=int(v[1])-1;end=start+len(v[3])
    if start>=b or end<=a:continue
    info=dict(item.split('=',1) for item in v[7].split(';') if '=' in item)
    key=(v[0],v[1],v[3],v[4]);found[key]=(l,dict(genes=w['genes'],chrom=chrom,start=start,end=end,ref=v[3],alt=v[4],filter=v[6],AF=info.get('AF'),AC=info.get('AC'),AN=info.get('AN'),qual=v[5]))
  raw.extend(x[0] for x in found.values());out.extend(x[1] for x in found.values())
  summary=dict(callset=name,genes=w['genes'],sample_count=nsamples,records=len(found),PASS=sum(x[1]['filter']=='PASS' for x in found.values()));allsummary.append(summary);print(summary,flush=True)
 if out:
  with (O/(name+'.local.tsv')).open('w',newline='') as f:
   wr=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');wr.writeheader();wr.writerows(out)
 (O/(name+'.local.vcf')).write_text('\n'.join(headers+raw)+'\n')
 (O/(name+'.provenance.json')).write_text(json.dumps(dict(url=base+name,index_sha256=hashlib.sha256(idx).hexdigest(),requests=requests,windows=windows),indent=2))
(O/'small_variant_summary.json').write_text(json.dumps(allsummary,indent=2))
