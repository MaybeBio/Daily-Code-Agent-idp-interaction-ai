"""Download every file listed in metadata/sources.tsv and check/record its sha256.

Run from the project folder:  pixi run download
"""
import csv
import hashlib
import urllib.request
from pathlib import Path

SOURCES = Path("metadata/sources.tsv")
CHECKSUMS = Path("metadata/checksums.sha256")

# checksums from earlier runs (file -> sha256)
known = {}
if CHECKSUMS.exists():
    for line in CHECKSUMS.read_text().splitlines():
        if line.strip():
            sha, name = line.split()
            known[name] = sha

with SOURCES.open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        path = Path(row["file"])
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.exists():
            print(f"downloading {path}")
            urllib.request.urlretrieve(row["url"], path)
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        name = path.as_posix()
        if name in known and known[name] != sha:
            raise SystemExit(f"checksum mismatch: {name} (file changed since it was recorded)")
        known[name] = sha

CHECKSUMS.write_text("".join(f"{sha}  {name}\n" for name, sha in sorted(known.items())))
print(f"{len(known)} files OK, checksums in {CHECKSUMS}")
