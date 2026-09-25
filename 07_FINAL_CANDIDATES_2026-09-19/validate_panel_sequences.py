from pathlib import Path
import pysam,csv,json,hashlib
P=Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino');O=P/'07_FINAL_CANDIDATES_2026-09-19'
fa=pysam.FastaFile(str(P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna'));genes=[]
for l in (P/'05_SHIP/canine_all_genes_stranded.bed').read_text().splitlines():
 v=l.split('\t')
 if len(v)>=6:genes.append(v)
rows=list(csv.DictReader((O/'shortlist_evidence.tsv').open(),delimiter='\t'))+list(csv.DictReader((O/'reserve_evidence.tsv').open(),delimiter='\t'));summary=[];core=[];flanks=[]
for r in rows:
 c=r['chrom'];s=int(r['start']);e=int(r['end']);seq=fa.fetch(c,s,e).upper();assert len(seq)==1000 and set(seq)<=set('ACGT')
 local=[v for v in genes if v[0]==c];dist=lambda v:max(int(v[1])-e,s-int(v[2]),0);d=min(map(dist,local));nearest=sorted({v[3] for v in local if dist(v)==d});assert d==int(r['gene_body_gap_bp']) and r['nearest_gene_body'] in nearest
 overlap=sum(max(0,min(e,int(v[2]))-max(s,int(v[1]))) for v in local);assert overlap==0
 flank=fa.fetch(c,s-1000,e+1000).upper();assert len(flank)==3000 and flank[1000:2000]==seq
 core.append(f">{r['window_id']} ROS_Cfam_1.0 {c}:{s}-{e} 0based_halfopen reference_forward decision={r['decision']}\n"+'\n'.join(seq[i:i+80] for i in range(0,1000,80))+'\n')
 flanks.append(f">{r['window_id']} ROS_Cfam_1.0 {c}:{s-1000}-{e+1000} 0based_halfopen reference_forward CONTEXT_NOT_VALIDATED_TARGET\n"+'\n'.join(flank[i:i+80] for i in range(0,3000,80))+'\n')
 summary.append(dict(window_id=r['window_id'],length=len(seq),N_count=seq.count('N'),GC_fraction=(seq.count('G')+seq.count('C'))/len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),reference_gene_body_overlap_bp=overlap,independent_nearest_gene_body=nearest,independent_gene_body_gap_bp=d,flank_length=3000,flank_sha256=hashlib.sha256(flank.encode()).hexdigest(),scope='Reference sequence only; no donor genotype, no guide sequence and no cut position selected'))
(O/'panel_plus_reserve_ROS_reference.fa').write_text(''.join(core));(O/'panel_plus_reserve_CONTEXT_flanks.fa').write_text(''.join(flanks));(O/'sequence_validation.json').write_text(json.dumps(dict(status='passed',windows=summary,reference='GCF_014441545.1 ROS_Cfam_1.0',flanks='Context only, not screened as safe intervals'),indent=2));print(json.dumps(summary,indent=2))
