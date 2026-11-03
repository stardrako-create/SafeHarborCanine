from pathlib import Path
import subprocess,json,hashlib,datetime,pysam
N=Path(__file__).resolve().parent;S=N/'ANO2_two_guides_full_UU_20261009';assert json.loads((S/'status.json').read_text())['state']=='completed'
I=Path('/mnt/d/Jin2024_work/bowtie1_UU_20261009');O=N/'full_UU_index_integrity_20261009_v2';O.mkdir(exist_ok=False);cmd=['/mnt/d/Jin2024_work/bowtie1_probe_20261009/extracted/usr/bin/bowtie-inspect','-s',str(I/'UU')];out=subprocess.check_output(cmd,text=True);(O/'index.summary.txt').write_text(out);(O/'command.json').write_text(json.dumps(cmd));names={}
for line in out.splitlines():
 x=line.split('\t')
 if x[0].startswith('Sequence-'):
  assert len(x)==3; ident=x[1].split()[0]; assert ident not in names; names[ident]=int(x[2])
fa=pysam.FastaFile('/mnt/d/Jin2024_work/uu_gsd/genome_uu_gsd.fa');expected=dict(zip(fa.references,fa.lengths));assert names==expected,'Index reference names/lengths mismatch'
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
files=list(I.glob('UU*.ebwt'));assert len(files)==6 and all(p.stat().st_size>0 for p in files)
result=dict(state='completed',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),references=len(names),total_reference_bp=sum(names.values()),index_matches_all_FASTA_reference_names_and_lengths=True,files={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files},summary_sha256=sha(O/'index.summary.txt'),script_sha256=sha(Path(__file__)),FASTA_provenance_receipt_sha256=sha(N/'full_UU_sequence_provenance_20261009.json'),limits=['Index structural/provenance verification plus alignment sequence checks; no independent biological confirmation.'])
(O/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k] for k in ['state','references','total_reference_bp']}))

