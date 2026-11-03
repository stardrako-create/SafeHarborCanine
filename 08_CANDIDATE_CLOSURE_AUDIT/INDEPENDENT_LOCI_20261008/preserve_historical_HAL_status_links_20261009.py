from pathlib import Path
N=Path(__file__).resolve().parent;P=N.parent.parent;files=[P/'CURRENT.md',N.parent/'CURRENT_STATE_AND_DECISIONS.md',Path('/mnt/d/Biblioteca/Biblioteca/Labs/Vasco Barreto Lab/Notas/Safe Harbor CAR-T - ESTADO ATUAL E DECISOES.md')];latest='## 09/10/2026 — retoma após interrupção e painéis de especificidade aprofundados';separator='### Histórico preservado a partir daqui'
for p in files:
 t=p.read_text();assert latest in t and separator in t;front,history=t.split(separator,1);prior,newnote=history.rsplit(latest,1);prior=prior.replace('remote_HAL_verified_20261009_v3/status.json','remote_HAL_verified_20261009_v2/status.json');p.write_text(front+separator+prior+latest+newnote)
print('Historical v2 provenance restored; current front and latest v3 note retained.')
