from pathlib import Path
import runpy, json, time, traceback
p=Path(__file__).parent
s=p/'SV_PBGV000010/alignment_chunked/continuation_status.json'
def status(stage, **kw):
 s.write_text(json.dumps(dict(stage=stage,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**kw),indent=2))
try:
 status('alignment_recovery_running')
 if not (s.parent/'completed.json').exists():
  runpy.run_path(str(p/'recover_sv_chunks.py'),run_name='__main__')
 status('descriptive_coverage_running')
 runpy.run_path(str(p/'summarize_sv_chunked_alignment.py'),run_name='__main__')
 status('descriptive_coverage_complete_pending_scientific_review')
except BaseException as e:
 status('failed_preserved',error=repr(e));raise
