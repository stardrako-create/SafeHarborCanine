from pathlib import Path
N=Path(__file__).resolve().parent
s=(N/'extract_verified_HAL_v2.py').read_text();s=s.replace('import http.server','import requests\nimport http.server') if 'import http.server' in s else 'import requests\n'+s
s=s.replace('remote_HAL_verified_20261009_v2','remote_HAL_verified_20261009_v3').replace('safeharbor_hal_verified_20261009_v2','safeharbor_hal_verified_20261009_v3').replace('268435456','536870912')
needle="""     req=urllib.request.Request(URL,headers={'Range':f'bytes={begin}-{end}'})
     with urllib.request.urlopen(req,timeout=30) as resp:
      if resp.status!=206 or resp.headers.get('Content-Range')!=f'bytes {begin}-{end}/{total}':raise RuntimeError('Invalid upstream range')
      assert resp.headers.get('ETag')==integrity['source']['ETag'];data=resp.read(end-begin+2)"""
replacement="""     with upstream_session.get(URL,headers={'Range':f'bytes={begin}-{end}'},stream=True,timeout=30) as resp:
      if resp.status_code!=206 or resp.headers.get('Content-Range')!=f'bytes {begin}-{end}/{total}':raise RuntimeError('Invalid upstream range')
      assert resp.headers.get('ETag')==integrity['source']['ETag'];data=resp.raw.read(end-begin+2)"""
assert needle in s;s=s.replace(needle,replacement);s=s.replace('class Handler(http.server.BaseHTTPRequestHandler):','upstream_session=requests.Session()\nclass Handler(http.server.BaseHTTPRequestHandler):');out=N/'extract_verified_HAL_v3.py';assert not out.exists();out.write_text(s)
s=(N/'score_verified_MAF_v2.py').read_text().replace('remote_HAL_verified_20261009_v2','remote_HAL_verified_20261009_v3').replace('remote_conservation_verified_20261009_v2','remote_conservation_verified_20261009_v3');out=N/'score_verified_MAF_v3.py';assert not out.exists();out.write_text(s);print('Prepared v3 only; no worker started')
