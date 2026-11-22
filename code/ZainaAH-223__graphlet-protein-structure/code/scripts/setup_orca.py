"""Download ORCA source at a pinned commit, check checksums, compile, self-test.

Run from the project folder:  pixi run setup-orca
"""
import filecmp
import hashlib
import os
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

COMMIT = "e146a8b1a99a90f5e3096b7bcc2ab0ea246c3ca7"
BASE = f"https://raw.githubusercontent.com/thocevar/orca/{COMMIT}/"
FILES = {  # file -> expected sha256
    "orca.cpp": "01c580febfd3653a632a04d6b466f7c6f2bf2f9baf3a75a2e553be30ac7ee778",
    "graph.in": "8d541fe333dc4fefa7993e429a47ecb7b2da0e6ce45b637702e76998b7aaadd6",
    "node-4.out": "0b1a605505ee3a5e974a50fff0da60f6aa87975dc779caddf3f61ec6b9e8743c",
}
OUT = Path("tools/orca")
OUT.mkdir(parents=True, exist_ok=True)

for name, sha in FILES.items():
    path = OUT / name
    if not path.exists():
        urllib.request.urlretrieve(BASE + name, path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != sha:
        sys.exit(f"checksum mismatch for {name}")

# CXX is set by Pixi's cxx-compiler package; fall back to system compilers
cxx = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
binary = OUT / "orca"
subprocess.run([cxx, "-O2", "-std=c++11", "-o", str(binary), str(OUT / "orca.cpp")], check=True)

# self-test: ORCA's bundled example must give ORCA's bundled answer
subprocess.run([str(binary), "node", "4", str(OUT / "graph.in"), str(OUT / "test.out")],
               check=True, capture_output=True)
if not filecmp.cmp(OUT / "test.out", OUT / "node-4.out", shallow=False):
    sys.exit("ORCA self-test failed")
print(f"ORCA OK (commit {COMMIT[:7]}, compiler {Path(cxx).name})")
