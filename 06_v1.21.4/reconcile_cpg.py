from pathlib import Path
import shutil,hashlib,json
P=Path(r'D:\Biblioteca\Bioquímica 1º ano 2025-2026\2º Semestre\Vasco M. Barreto\Safe Harbor Canino')
O=P/'06_v1.21.4';(O/'code').mkdir(exist_ok=True);(O/'prior_main').mkdir(exist_ok=True)
p=P/'scripts/build_cpg_islands.py';text=p.read_text(encoding='utf-8')
shutil.copy2(p,O/'prior_main/build_cpg_islands.py')
start=text.index('        # A region called by both parameterizations')
end=text.index('        print(f"{chrom}:',start)
text=text[:start]+'''        # Preserve each definition's own coordinates and statistics.
        rows.extend(independent_calls(chrom, per_type_regions))

'''+text[end:]
marker='def main():'
assert marker in text
text=text.replace(marker,'''def independent_calls(chrom, per_type_regions):
    """Overlap alone never implies that an interval satisfies both definitions."""
    return [(chrom, s, e, name, length, gc, oe, n)
            for name, regions in per_type_regions.items()
            for s, e, length, gc, oe, n in regions]


def main():''')
text=text.replace('A region satisfying both is emitted once, typed CGI_GGF+CGI_TJ.','Each definition is emitted independently with its own coordinates and statistics.\nOverlapping GGF and TJ calls are retained; no overlap-derived dual labels.')
text=text.replace('    with open(args.out_bed, "w", encoding="utf-8") as f:', '    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))\n    with open(args.out_bed, "w", encoding="utf-8") as f:')
text=text.replace('{n_ggf} GGF-only, {n_tj} TJ-only','{n_ggf} GGF calls, {n_tj} TJ calls (independent)')
p.write_text(text,encoding='utf-8');shutil.copy2(p,O/'code/build_cpg_islands.py')
old=P/'05_SHIP/cpg_islands_ROS_Cfam_1.0.bed';shutil.copy2(old,O/'prior_main/cpg_islands_ROS_Cfam_1.0.bed')
corrected=P/'06_v1.21.1/cpg_islands_independent_calls.bed'
shutil.copy2(corrected,old)
(O/'cpg_reconciliation.json').write_text(json.dumps({'corrected_main_script':str(p),'main_bed':str(old),'reused_validated_genomewide_calls':str(corrected),'main_bed_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),'source_bed_sha256':hashlib.sha256(corrected.read_bytes()).hexdigest(),'historical_backup':str(O/'prior_main'),'new_full_genome_build':False},indent=2),encoding='utf-8')
print('Main script corrected; validated BED promoted with historical backup')
