"""build_original.py: copies every script of audit/manifest.csv into original/, byte for byte, keeping each project's folder layout,
and writes original/SHA256SUMS so that `cd original && shasum -a 256 -c SHA256SUMS` proves the copies are intact. The BENDER-MT
files are left out on purpose; helpers/ holds the definitions the paper needs (tools/extract_helpers.py).
    python3 -I tools/build_original.py
"""
import csv
import hashlib
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "original"
PREFIX = {"/Volumes/User Homes/Nu Project/BENDER-DYN/redo_evidence/": "redo_evidence/",
          "/Volumes/User Homes/Nu Project/BENDER-DYN/methods_paper/figures/": "methods_paper/figures/",
          str(Path.home() / "Downloads/Scripts/bender_reanalysis") + "/": "bender_reanalysis/"}
SKIP = "/Volumes/User Homes/Nu Project/BENDER-MT/"


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    sums, skipped, unknown = [], [], []
    for r in csv.DictReader(open(ROOT / "audit/manifest.csv")):
        src = r["file"]
        if src.startswith(SKIP):
            skipped.append(Path(src).name)
            continue
        pre = next((p for p in PREFIX if src.startswith(p)), None)
        if pre is None:
            unknown.append(src)
            continue
        dst = OUT / PREFIX[pre] / src[len(pre):]
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        h = hashlib.sha256(dst.read_bytes()).hexdigest()
        assert h == r["sha256"], src
        sums.append(f"{h}  {dst.relative_to(OUT)}")
    (OUT / "SHA256SUMS").write_text("\n".join(sorted(sums, key=lambda s: s.split("  ")[1])) + "\n")
    print(f"copied {len(sums)} scripts; left out {skipped}; outside the known folders: {unknown}")


if __name__ == "__main__":
    main()
