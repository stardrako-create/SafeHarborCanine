#!/usr/bin/env python3
"""
Verify the files of a snapshot against its SHA-256 manifest (standard library only).

Usage (from the repository root or from an unpacked archive):
    python scripts/verify_manifest.py PUBLICATION_MANIFEST_20261012.json [--root .]

The manifest hashes the Git index bytes with LF line endings and excludes itself. Text files that were
checked out with CRLF are re-tried with the line endings normalised, so a Windows checkout verifies too.
Exit status is 0 when every listed file is present and matches, 1 otherwise. Files that exist in the tree
but are not listed (for example local-only material that .gitignore keeps out of the archive) are reported
as a count and do not fail the check.
"""
import argparse
import hashlib
import json
import os
import sys


def sha256_of(path, normalise_crlf=False):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        if normalise_crlf:
            h.update(f.read().replace(b'\r\n', b'\n'))
        else:
            for chunk in iter(lambda: f.read(1 << 20), b''):
                h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('manifest')
    ap.add_argument('--root', default='.')
    a = ap.parse_args()
    man = json.load(open(a.manifest, encoding='utf-8'))
    files = man['files']
    ok = missing = bad = 0
    problems = []
    for e in files:
        p = os.path.join(a.root, *e['path'].split('/'))
        if not os.path.isfile(p):
            missing += 1
            problems.append('MISSING  ' + e['path'])
            continue
        if sha256_of(p) == e['sha256'] or sha256_of(p, normalise_crlf=True) == e['sha256']:
            ok += 1
        else:
            bad += 1
            problems.append('MISMATCH ' + e['path'])
    listed = {e['path'] for e in files} | {os.path.basename(a.manifest)}
    extra = 0
    for dirpath, dirnames, filenames in os.walk(a.root):
        dirnames[:] = [d for d in dirnames if d != '.git']
        for fn in filenames:
            rel = os.path.relpath(os.path.join(dirpath, fn), a.root).replace(os.sep, '/')
            if rel not in listed:
                extra += 1
    print('snapshot %s: %d files listed | %d verified | %d missing | %d mismatching | %d files in the tree not listed' % (
        man.get('snapshot', '?'), len(files), ok, missing, bad, extra))
    for p in problems[:50]:
        print('  ' + p)
    if len(problems) > 50:
        print('  ... and %d more' % (len(problems) - 50))
    return 0 if not problems else 1


if __name__ == '__main__':
    sys.exit(main())
