from pathlib import Path
import pysam,json,hashlib
A=Path(__file__).resolve().parent;paths=['/mnt/d/SafeHarborWork/w11_spades_20261004_ascii/assembly_libfix/contigs.fasta','/mnt/d/SafeHarborWork/w11_augmented_20261005/assembly/contigs.fasta'];sets=[]
for p in paths:
 with pysam.FastxFile(p) as f:sets.append({r.sequence for r in f})
out={'identical_sequence_sets':sets[0]==sets[1],'old_only':len(sets[0]-sets[1]),'new_only':len(sets[1]-sets[0]),'sha256':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths}};(A/'W11_AUGMENTED_ASSEMBLY_20261005/comparison.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))