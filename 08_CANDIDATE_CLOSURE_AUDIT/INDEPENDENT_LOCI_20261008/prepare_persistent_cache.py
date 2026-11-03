from pathlib import Path
N=Path(__file__).resolve().parent
s=(N/'extract_remote_HAL_cached.py').read_text(encoding='utf-8-sig')
s=s.replace("remote_HAL_cached_20261009'","remote_HAL_persistent_20261009'").replace('safeharbor_hal_cached_20261009','safeharbor_hal_persistent_20261009')
s=s.replace("cache={};lock=threading.Lock();", "disk=Path('/mnt/d/Jin2024_work/hal_range_cache_20261009');disk.mkdir(exist_ok=True);cache={};lock=threading.Lock();")
s=s.replace("meta=dict(resp.headers)","meta=dict(resp.headers)\n    signature={k:meta.get(k) for k in ['Content-Length','Last-Modified','ETag']};signature['url']=URL\n    manifest=disk/'source.json'\n    if manifest.exists():assert json.loads(manifest.read_text())==signature,'Remote source changed'\n    else:manifest.write_text(json.dumps(signature,indent=2))")
s=s.replace("if begin not in cache:\n     with lock:","blockfile=disk/(str(begin)+'.bin')\n    if begin not in cache and blockfile.exists():\n     data=blockfile.read_bytes();assert len(data)==end-begin+1;cache[begin]=data\n    if begin not in cache:\n     with lock:")
s=s.replace("cache[begin]=data\n     with lock:","cache[begin]=data;temporary=blockfile.with_suffix('.tmp');temporary.write_bytes(data);temporary.replace(blockfile)\n     with lock:")
(N/'extract_remote_HAL_persistent.py').write_text(s)
print('Prepared only; not launched while existing extraction is active.')
