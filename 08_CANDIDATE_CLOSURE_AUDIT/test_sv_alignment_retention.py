import tempfile,subprocess,json
from pathlib import Path
S='/home/stardrako/miniforge3/envs/atac/bin/samtools'
with tempfile.TemporaryDirectory() as td:
 d=Path(td);sam=d/'fixture.sam';seq='A'*50;qual='I'*50
 lines=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:NC_049253.1\tLN:10000','@SQ\tSN:control\tLN:10000']
 for name,c1,c2 in [('target','NC_049253.1','NC_049253.1'),('cross','control','NC_049253.1'),('other','control','control')]:
  for flag,c,pos,mate,mpos in [(99,c1,101,c2,201),(147,c2,201,c1,101)]:
   lines.append(f'{name}\t{flag}\t{c}\t{pos}\t60\t50M\t{mate}\t{mpos}\t0\t{seq}\t{qual}')
 sam.write_text('\n'.join(lines)+'\n')
 def run(*args):return subprocess.run([S,*map(str,args)],check=True,capture_output=True,text=True).stdout
 run('view','-b','-e','rname == "NC_049253.1" || rnext == "NC_049253.1"','-o',d/'selected.bam',sam)
 run('sort','-n','-o',d/'name.bam',d/'selected.bam');run('fixmate','-m',d/'name.bam',d/'fix.bam');run('sort','-o',d/'coord.bam',d/'fix.bam');run('markdup',d/'coord.bam',d/'final.bam');run('quickcheck',d/'final.bam');run('index',d/'final.bam')
 records=run('view',d/'final.bam').splitlines();assert len(records)==4
 assert sorted(r.split('\t')[0] for r in records)==['cross','cross','target','target']
 print(json.dumps({'test':'regional retention includes both discordant mates and excludes unrelated pair','records':4,'sort_fixmate_markdup_index':'passed'}))
