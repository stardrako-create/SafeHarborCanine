from pathlib import Path
import json,hashlib,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent;A=N.parent/'ALL_CANDIDATES_20261005';old=A/'ANO2_REVIEW_AND_DEPENDENCIES_20261009.md';new=A/'ANO2_REVIEW_AND_DEPENDENCIES_20261009_v2.md';assert not new.exists()
s=old.read_text().replace('Um painel alternativo23nt+TTTV está em execução, como diagnóstico para apoiar a escolha.','O painel alternativo23nt+TTTV terminou para80sequências em3loci. EmANO2,5/14não têm hits no modelo canónico; ao alargar aPAMs com uma substituiçãoT>C,2/14continuam semhits no modelo definido. É um diagnóstico para apoiar a escolha.');s+='''

## Adendo: duas hipóteses de sequência com resultado computacional favorável

| Hipótese | Protospacer ROS BED0 | Strand | PAM | Resultado do modelo |
|---|---|---|---|---|
| ANO2_TTTV23_1 |NC_051831.1:39708151–39708174 |+ |TTTC |0hits fora do locus no modelo definido |
| ANO2_TTTV23_2 |NC_051831.1:39708436–39708459 |+ |TTTC |0hits fora do locus no modelo definido |

Modelo explícito:23nt; PAMsTTTV/CTTV/TCTV/TTCV; alinhamentosungapped até3substituições ou exatamenteumbulge de1base até2substituições; referênciaROS pesquisável. Foram testadas as duas orientações. Protospacer ePAM mantêm-se intactos nas32combinações sintéticas dos5indelsPASS locais. A construção e a pesquisa foram conferidas com exemplos de22/23/24bases contra enumeração exaustiva. Recibos e sequências em `ANO2_TTTV23_REFERENCE_MODEL_CANDIDATES_20261009.tsv/json` no diretório de análise.

**Isto não define a nuclease real nem aprova guias.** O modelo não abrange todas as especificidades dePAM de todas as variantes da enzima, múltiplosbulges, variantes emoutrosloci ou genomasindividuais. Basesambíguas/lacunas doassembly são informação emfalta, não prova deausência deoff-targets. Os alinhamentos comumbulge que sobrepõem o intervalo pretendido foram separados do rastreio de outros loci; não se demonstrou ausência decortes deslocados localmente. Zero neste modelo não significa zerooff-targets biológicos.

## Retoma após interrupção

09/10/2026: o utilizador informou que teve dedesligar. A tentativaHALv2 ficou comestado`running`,mas semprocessos presentes; a causa precisa no sistemaoperativo não foi estabelecida. OMAFparcial foi preservado semvalidação. Iniciado umúnico workerHALv3 destacado (launcherWindows27748,wrapperLinux292), usando cacheexistente e novosoutputs. Aumentado apenasocap operacional de transferência nova para512MiB; pedidosHTTPS mantêm validação206/Content-Range/ETag e gravação durável. Score só decorre após recibodeextração`completed`. Não alterados limiares científicos. Estado em `verified_HAL_v3_chain_20261009/status.json` e `remote_HAL_verified_20261009_v3/status.json`.
''';new.write_text(s)
files=[P/'CURRENT.md',N.parent/'CURRENT_STATE_AND_DECISIONS.md',N.parent/'CANDIDATE_EXPANSION_20261005.md',A/'CANDIDATO_ANO2_39708.md',N.parent/'MEMPALACE_SYNC_PENDING_20260922.md',Path('/mnt/d/Biblioteca/Biblioteca/Labs/Vasco Barreto Lab/Notas/Safe Harbor CAR-T - ESTADO ATUAL E DECISOES.md')]
note='''\n\n## 09/10/2026 — retoma após interrupção e painéis de especificidade aprofundados\nHALv2 estava semworker e comrecibo`running`desatualizado após desligar;MAFparcial preservado,nãovalidado. HALv3 único iniciado destacado, novosoutputs e cacheexistente; cap512MiB novos, mesmoscritérios. LauncherWindows27748/wrapperLinux292; verificarverified_HAL_v3_chain_20261009/status.json e remote_HAL_verified_20261009_v3/status.json antes derelançar.\nANO2:102protospacers+PAMNGG intactos em32combinaçõesindel;9perfisNGG revistos contraanotaçõesgénicas e nenhumaprovado. Painel23ntTTTV:14hipótesesANO2;5semhitsno modelocanónico;2mantêmzero aoalargar paraTTTV/CTTV/TCTV/TTCV,até3substituições ungapped ou1bulge até2substituições. Estes2alvos/PAM estãointactos em32combinações. Basesambíguas/genomasalternativos, atividadeenucleasereal nãocobertos;nãosão guiasebiologiavalidados. Nucleaserealcontinua perguntada, semresposta nestaexecução.\nRelatório novoANO2_REVIEW_AND_DEPENDENCIES_20261009_v2.md preservaov1eassuas fontes. Checkpoints anterioresnãoalterados;nãodeclararaprovaçãointegraldaauditoria.MemPalaceSSE404: filalocal,nãosync.\n'''
for p in files:
 with p.open('a') as f:f.write(note)
for p in [P/'CURRENT.md',N.parent/'CURRENT_STATE_AND_DECISIONS.md',files[-1]]:
 t=p.read_text();t=t.replace('extract_verified_HAL_v2.py (conservação3loci; scorer encadeado)','extract_verified_HAL_v3.py (conservação3loci; scorer encadeado via worker destacado)').replace('remote_HAL_verified_20261009_v2/status.json','remote_HAL_verified_20261009_v3/status.json');p.write_text(t)
sourcepaths=[new,N/'hypothetical_TTTV23_alternate_PAM_sensitivity_20261009.json',N/'ANO2_TTTV23_REFERENCE_MODEL_CANDIDATES_20261009.json',N/'hypothetical_TTTV23_search_validation_20261009.json',N/'ANO2_hypothetical_TTTV23_PAM_indel_robustness_20261009.json',N/'HAL_interruption_recovery_20261009.json']
(N/'ANO2_review_delivery_manifest_20261009_v2.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),complete_audit_approval=False,biological_validation=False,files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in sourcepaths]),indent=2));print(str(new))
