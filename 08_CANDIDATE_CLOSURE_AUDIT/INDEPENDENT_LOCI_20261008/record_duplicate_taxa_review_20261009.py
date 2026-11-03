from pathlib import Path
import json,hashlib,datetime
N=Path(__file__).resolve().parent;P=N.parent.parent;A=N.parent/'ALL_CANDIDATES_20261005';receipt=N/'ANO2_duplicate_taxa_sensitivity_20261009_v2/receipt.json';x=json.loads(receipt.read_text());assert x['state']=='completed' and all(not r['exceeds_legacy_6_5'] for r in x['results']) and x['masked_nonreference_ACGT']['zero_refbases']==0
failed=N/'conservation_alignment_information_review_20261009_failed.json';failed.write_text(json.dumps(dict(state='failed',source_script='review_conservation_alignment_information.py',reason='Duplicate species-prefix rows in retained ANO2MAF; initial assertion stopped before resultapproval',superseded_by=str(receipt),original_inputs_preserved=True),indent=2))
old=A/'ANO2_REVIEW_AND_DEPENDENCIES_20261009_v2.md';new=A/'ANO2_REVIEW_AND_DEPENDENCIES_20261009_v3.md';assert not new.exists();s=old.read_text();s+='''

## Qualidade do alinhamento: duplicações de espécies resolvidas na sensibilidade

A verificação adicional encontrou espécies com vários scaffolds no mesmo bloco do MAFretido. Em4105pares espécie/coluna de referência comduplicação,47contêm basesACGT discordantes. Inverter a ordem dascópias altera65scores porbase: a escolha dolinha importa no processamento instalado. Não tratar estas cópias como confirmação independente.

Reanálise conservadora, predefinida seminspecionar scores para escolherlinhas: retirar todas aslinhas dequalquer espécie não referência que esteja duplicada dentro do bloco. DOG e coordenadas permanecem intactos. Os1000scoresROS mantêm cobertura; cada base tem2–194espécies não referência comACGT, mediana164. Não há baseDOG seminformaçãoACGT de outraespécie neste alinhamento curado. Isto descreve suporte desigual, não impõe um novo limiar deespécies.

Máximo da média50bp: original1,56426; ordem de cópias invertida1,54716; espéciesduplicadas retiradas1,55246. Todas asanálises ficam abaixo dolimiar anterior6,5. A conclusão do gate é robusta à ambiguidade examinada; o valor curado1,55246 deve acompanhar as limitações. Não mudar valores históricos nem chamar a estasreanálises dadosindependentes.

Recibo: `../INDEPENDENT_LOCI_20261008/ANO2_duplicate_taxa_sensitivity_20261009_v2/receipt.json`; MAF/WIG/comandos/hashes separados. A tentativa deQA inicial parou naassertiondeduplicações e está registada em `conservation_alignment_information_review_20261009_failed.json`; não emitiu resultado aprovado. A tentativa de sensibilidade inicial parou ao reconhecer a linhaMAF`a`; corrigida noutro script/diretório, mantendo osanteriores.
''';new.write_text(s)
note='''\n\n## 09/10/2026 — controlo de conservação com espécies duplicadas\nMAFretidoANO2:4105pares espécie/coluna duplicados,47comACGT discordante. Inverter ordem altera65scores; retirartodaslinhasdeespéciesduplicadasporbloco dámax50=1,55246, comparado comraw1,56426 eordeminvertida1,54716;todos<6,5. Nãoescolhidaslinhasporscore.1000basesROScuradas,min2/mediana164/max194espéciesnãoDOGcomACGT;zerobasesseminformaçãodeoutraespécie. ReciboANO2_duplicate_taxa_sensitivity_20261009_v2/receipt.json. QAinicialfalhouporduplicações,registado,semaprovação;inputseoutputshistóricospreservados. RelatórionovoANO2_REVIEW_AND_DEPENDENCIES_20261009_v3.md. Nãoévalidaçãobiológica nemaprovaçãointegraldaintervenção.\n'''
files=[P/'CURRENT.md',N.parent/'CURRENT_STATE_AND_DECISIONS.md',N.parent/'CANDIDATE_EXPANSION_20261005.md',A/'CANDIDATO_ANO2_39708.md',N.parent/'MEMPALACE_SYNC_PENDING_20260922.md',Path('/mnt/d/Biblioteca/Biblioteca/Labs/Vasco Barreto Lab/Notas/Safe Harbor CAR-T - ESTADO ATUAL E DECISOES.md')]
for p in files:
 with p.open('a') as f:f.write(note)
for p in [P/'CURRENT.md',N.parent/'CURRENT_STATE_AND_DECISIONS.md',files[-1]]:
 t=p.read_text();front,history=t.split('### Histórico preservado a partir daqui',1);front=front.replace('Conservação1000bp recalculada,max50=1,56426<6,5','Conservação1000bp curada,max50=1,55246<6,5, robusta à ordem/remoção de espécies duplicadas');p.write_text(front+'### Histórico preservado a partir daqui'+history)
(N/'ANO2_review_delivery_manifest_20261009_v3.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),complete_audit_approval=False,biological_validation=False,files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in [new,receipt,failed]]),indent=2));print(str(new))
