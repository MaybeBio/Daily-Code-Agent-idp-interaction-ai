"""build_original.py: copies every included file of audit/manifest.csv into original/, byte for byte, and writes original/SHA256SUMS,
so that `cd original && shasum -a 256 -c SHA256SUMS` proves the copies are intact.
    python3 -I tools/build_original.py
"""
import csv
import hashlib
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "original"
PREFIX = {"/Volumes/User Homes/Nu Project/BENDER-MT/code/": "code/",
          str(Path.home() / "Downloads/Scripts/bender_reanalysis") + "/": "bender_reanalysis/"}


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    sums = []
    for r in csv.DictReader(open(ROOT / "audit/manifest.csv")):
        if r["included"] != "yes":
            continue
        pre = next(p for p in PREFIX if r["file"].startswith(p))
        dst = OUT / PREFIX[pre] / r["file"][len(pre):]
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(r["file"], dst)
        h = hashlib.sha256(dst.read_bytes()).hexdigest()
        assert h == r["sha256"], r["file"]
        sums.append(f"{h}  {dst.relative_to(OUT)}")
    (OUT / "SHA256SUMS").write_text("\n".join(sorted(sums, key=lambda s: s.split("  ")[1])) + "\n")
    print(f"copied {len(sums)} files")


if __name__ == "__main__":
    main()
