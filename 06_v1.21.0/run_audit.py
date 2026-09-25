"""Isolated v1.21.0 scoring replay and evidence audit; run under the existing WSL atac environment."""
from pathlib import Path
import csv, hashlib, json, shutil, subprocess, sys, datetime, collections
import pyBigWig

P = Path('/mnt/d/Biblioteca/Bioquímica 1º ano 2025-2026/2º Semestre/Vasco M. Barreto/Safe Harbor Canino')
O = P / '06_v1.21.0'

def read(path):
    with path.open(encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def key(r):
    return r['chrom'], int(r['start']), int(r['end'])

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024), b''): h.update(block)
    return h.hexdigest()

def write_table(path, rows):
    with path.open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)

def classify(hard_veto, missing):
    return 'excluded' if hard_veto else ('insufficient_evidence' if missing else 'passes_recorded_checks')

def main():
    # Regression assertions: unknown evidence must never be promoted to pass.
    assert classify(False, ['repeat']) == 'insufficient_evidence'
    assert classify(True, ['repeat']) == 'excluded'
    assert classify(False, []) == 'passes_recorded_checks'
    if O.exists(): raise SystemExit('Output exists: refusing to overwrite a run')
    O.mkdir(); (O/'snapshot').mkdir()
    for name in ['score_ship_candidates_v2.py','bw_utils.py','extract_genes_stranded.py','build_cpg_islands.py','sensitivity_analysis.py']:
        shutil.copy2(P/'scripts'/name, O/'snapshot'/name)
    shutil.copy2(Path(__file__), O/'run_audit.py')
    sources = {
        'candidates': '05_SHIP/ship_raw_candidates.tsv',
        'atac-mean-bw':'06_v1.20.0/ATAC/full76/mother_track_accessibility_level.bw',
        'atac-variability-bw':'06_v1.20.0/ATAC/full76/variability.bw',
        'atac-peaks-bed':'04_tracks_processadas/ROS_Cfam_1.0/ATAC/mother_track/consensus_peaks_ATAC.bed',
        'atac-peak-frequency-bw':'06_v1.20.0/ATAC/full76/peak_frequency.bw',
        'rrbs-mean-bw':'06_v1.20.0/RRBS/methylation_weighted_mean.bw',
        'rrbs-variability-bw':'06_v1.20.0/RRBS/variability.bw',
        'rrbs-coverage-bw':'06_v1.20.0/RRBS/cpg_coverage_frequency.bw',
        'tad-boundaries-bed':'04_tracks_processadas/ROS_Cfam_1.0/HiC/tad_boundaries.bed',
        'risk-genes':'05_SHIP/canine_risk_genes.tsv',
        'mappability-tsv':'05_SHIP/mappability_check.tsv',
        'mirna-bed':'05_SHIP/canine_miRNA.bed',
        'all-genes-bed':'05_SHIP/canine_all_genes.bed',
        'all-genes-stranded-bed':'05_SHIP/canine_all_genes_stranded.bed',
        'lncrna-smallrna-bed':'05_SHIP/canine_lncRNA_smallRNA.bed',
        'tad-intervals-bed':'05_SHIP/canine_tad_intervals.bed',
        'repeat-content-tsv':'05_SHIP/repeat_content_v5candidates.tsv',
        'ultraconserved-tsv':'05_SHIP/phylop_ultraconserved_v3_pilot_forscoring.tsv',
        'external-regulatory-bed':'05_SHIP/ehsan_regulatory_elements_ROS.bed',
    }
    candidates=read(P/sources['candidates'])
    assert len(candidates)==461 and len({key(r) for r in candidates})==461
    fai=P/'01_referencia/ROS_Cfam_1.0/GCF_014441545.1_ROS_Cfam_1.0_genomic.fna.fai'
    sizes={f[0]:int(f[1]) for line in fai.read_text().splitlines() if (f:=line.split('\t'))}
    # Validate all requested ranges against the reference and each input track.
    for arg,rel in sources.items():
        if not (P/rel).exists(): raise FileNotFoundError(rel)
        if rel.endswith('.bw'):
            with pyBigWig.open(str(P/rel)) as bw:
                for chrom,s,e in map(key,candidates):
                    assert bw.chroms(chrom)==sizes[chrom] and 0<=s<e<=sizes[chrom], (arg,chrom)
    manifest={'version':'1.21.0','scope':'scoring replay of existing tracks plus explicit evidence states; not a raw-data rebuild',
              'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'inputs':{},'source_snapshot':{}}
    for arg,rel in sources.items():
        f=P/rel
        manifest['inputs'][arg]={'path':str(f),'size':f.stat().st_size,'mtime_ns':f.stat().st_mtime_ns,'sha256':digest(f)}
    for f in (O/'snapshot').iterdir(): manifest['source_snapshot'][f.name]=digest(f)
    command=[sys.executable,str(O/'snapshot/score_ship_candidates_v2.py')]
    for arg,rel in sources.items(): command += ['--'+arg,str(P/rel)]
    command += ['--min-atac-accessibility-percentile','0.55','--out-scored',str(O/'raw_scoring.tsv'),'--out-passing-bed',str(O/'provisional_passing.bed')]
    manifest['command']=command
    (O/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print('Running frozen scorer with validated inputs',flush=True)
    with (O/'scoring.log').open('w') as log: subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    rows=read(O/'raw_scoring.tsv'); old={key(r):r for r in read(P/'06_v1.20.0/candidates_scored_v120d.tsv')}
    differences=[]
    for r in rows:
        for k in ['hard_veto','final_score','pct_repeat','max_50bp_rolling_phyloP','gene_5prime_clearance']:
            if r[k]!=old[key(r)][k]: differences.append((key(r),k,old[key(r)][k],r[k]))
    (O/'replay_comparison.json').write_text(json.dumps(differences,indent=2),encoding='utf-8')
    maps={key(r) for r in read(P/sources['mappability-tsv'])}
    tads=collections.defaultdict(list)
    for line in (P/sources['tad-intervals-bed']).read_text().splitlines():
        c,s,e=line.split('\t')[:3];tads[c].append((int(s),int(e)))
    diagnostics=[]
    handles={arg:pyBigWig.open(str(P/rel)) for arg,rel in sources.items() if rel.endswith('.bw')}
    for r in rows:
        c,s,e=key(r);missing=[]
        for field in ['pct_repeat','max_50bp_rolling_phyloP']:
            if not r[field]: missing.append(field)
        if key(r) not in maps: missing.append('mappability')
        if not any(a <= (s+e)//2 < b for a,b in tads[c]):missing.append('TAD_assignment')
        if r['no_rrbs_coverage']=='True':missing.append('RRBS')
        r['evidence_missing']=';'.join(missing)
        r['audit_status']=classify(r['hard_veto']=='True',missing)
        if r['hard_veto']=='False':
            d={'chrom':c,'start':s,'end':e,'genes':r['left_gene']+'/'+r['right_gene'],'score':r['final_score'],'status':r['audit_status'],'missing':r['evidence_missing']}
            for arg,bw in handles.items():
                d[arg+'_max']=bw.stats(c,s,e,type='max',exact=True)[0]
                d[arg+'_covered_fraction']=bw.stats(c,s,e,type='coverage',exact=True)[0]
            diagnostics.append(d)
    for bw in handles.values():bw.close()
    write_table(O/'candidates_scored_v121.tsv',rows);write_table(O/'provisional_candidate_diagnostics.tsv',diagnostics)
    with (O/'candidates_passing_recorded_checks.bed').open('w') as f:
        for r in sorted(rows,key=lambda r:float(r['final_score']),reverse=True):
            if r['audit_status']=='passes_recorded_checks':f.write('\t'.join(map(str,(*key(r),r['left_gene']+'/'+r['right_gene'],r['final_score'])))+'\n')
    cpg=read(P/'05_SHIP/cpg_islands_ROS_Cfam_1.0.bed')
    bad=[r for r in cpg if r['type']=='CGI_GGF+CGI_TJ' and (int(r['length'])<500 or float(r['gc_frac'])<.55 or float(r['obs_exp'])<.65)]
    write_table(O/'cpg_dual_label_threshold_failures.tsv',bad)
    counts=collections.Counter(r['audit_status'] for r in rows)
    summary={'n':len(rows),'status_counts':dict(counts),'v120d_replay_differences':len(differences),'cpg_dual_label_threshold_failures':len(bad),'diagnostics':diagnostics}
    (O/'audit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    report=f'''# Auditoria independente — v1.21.0, 2026-09-12

Execução do scoring sobre as tracks existentes v1.20.0, mantendo p55, vetos e pesos. Esta versão acrescenta classificação explícita de evidência ausente, verificações de assembly/intervalos, hashes e diagnósticos locais. Não reconstruiu tracks nem gerou novos dados de conservação.

- Candidatos: {len(rows)}.
- Estados de evidência: {dict(counts)}.
- Diferenças com v120d nos campos comparados: {len(differences)} (replay_comparison.json).
- Linhas CpG com etiqueta dupla que não satisfazem TJ nas coordenadas reportadas: {len(bad)}.

## Resultados materiais

1. Os três novos sobreviventes da v120d têm repetições e conservação por avaliar. A classificação anterior promovia ausência de dados a ausência de risco. Aqui ficam como insufficient_evidence, não como candidatos reprovados biologicamente.
2. Etiqueta CGI_GGF+CGI_TJ significa no código apenas sobreposição entre duas chamadas, não que o intervalo exportado satisfaça ambos os critérios. Preservar os dois intervalos e os seus identificadores/estatísticas é necessário.
3. 55% de intervalos externos a sobrepor CpG não se compara diretamente com a fração de bases do genoma em CpG. É necessário um modelo nulo ao nível dos intervalos; sobreposição não prova atividade de promotor.
4. As distâncias de 50 kb usam extremidades da janela; as de 300 kb usam o centro. São geometrias diferentes. A regra de 50 kb não define o futuro ponto de inserção. Atribuições a TADs são aproximações entre boundaries, não validação funcional.
5. O snapshot do builder não foi aplicado às tracks: estas ainda antecedem a recuperação de blocos corrompidos. Os diagnósticos aqui apresentados refletem as tracks existentes.
6. O relatório de sensibilidade antigo refere-se à v120c, não aos novos quatro da v120d. Percentagens de configurações não são probabilidades de segurança ou sucesso.

## Limites

passes_recorded_checks não significa safe harbor validado nem validação dos inputs de conservação/risco. Máscaras de cobertura não medem por si suporte biológico suficiente entre cães. CpG continua anotação complementar, sem novo veto. Não foram enviados emails, alteradas versões anteriores ou publicadas alterações.

## Reprodução

manifest.json contém o comando, inputs, hashes e runtime. scoring.log contém a saída completa. raw_scoring.tsv conserva a saída original; candidates_scored_v121.tsv acrescenta evidence_missing e audit_status. Os ficheiros foram separados para permitir auditar a mudança de interpretação sem alterar scores silenciosamente.
'''
    (O/'AUDIT.md').write_text(report,encoding='utf-8')
    (O/'CURRENT.md').write_text('# v1.21.0 — experimental audit run\n\nCompleted scoring replay; provisional candidates require missing evidence. Read AUDIT.md and manifest.json. Not a production replacement or full track rebuild.\n',encoding='utf-8')
    # Detect changes during the run rather than silently attributing results to a moving input.
    for arg,rel in sources.items():
        f=P/rel
        assert f.stat().st_mtime_ns==manifest['inputs'][arg]['mtime_ns'] and f.stat().st_size==manifest['inputs'][arg]['size'], 'Input changed: '+rel
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
