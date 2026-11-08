"""verify_frozen.py — check every entry of FROZEN.sha256 against the file on disk. Entries name files relative to this
folder (gaivi/); the two oldest model entries were written relative to redo_evidence/, so both bases are tried. An entry
whose file was later amended and re-hashed (a later entry for the same file matches) is reported as superseded, which is
not a problem; any other mismatch is.
  ~/.venvs/bender/bin/python verify_frozen.py"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
entries = []
for line in open(os.path.join(HERE, 'FROZEN.sha256')):
    if line.strip():
        digest, name = line.split()[:2]
        path = next((os.path.normpath(os.path.join(b, name)) for b in (HERE, os.path.dirname(HERE))
                     if os.path.isfile(os.path.join(b, name))), None)
        entries.append((digest, name, path))
actual = {p: hashlib.sha256(open(p, 'rb').read()).hexdigest() for _, _, p in entries if p}
bad = 0
for k, (digest, name, path) in enumerate(entries):
    if path is None:
        status = 'MISSING'
    elif actual[path] == digest:
        status = 'OK'
    elif any(p == path and d == actual[path] for d, _, p in entries[k + 1:]):
        status = 'SUPERSEDED (amended and re-hashed later)'
    else:
        status = 'CHANGED'
    bad += status in ('MISSING', 'CHANGED')
    print(f'{status:42s} {name}')
print(f'{bad} problem(s)')
sys.exit(1 if bad else 0)
