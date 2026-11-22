#!/usr/bin/env python3
"""
Step 11B (Final Phase 1 Exact Parity Lock):
  1. Adds IDR-immune NIS guard `if (hasOrdCore && platConf(r) < 0.35) continue;` to index.html.
  2. Adds SS (Disulfide) and Clash tags to renderContactInspectorTable() in index.html.
  3. Removes legacy `if min(len(residues_a), len(residues_b)) >= 50:` guard in unippi_cli.py.
  4. Validates syntax via ast.parse() and re-runs Test A (fold_vip1_bzip61.zip) & Test B (--generate-demo).
"""

import ast
import subprocess
import sys
from pathlib import Path

idx_path = Path("index.html")
cli_path = Path("unippi_cli.py")

# 1. Patch index.html: IDR-immune NIS + SS/Clash tags in Contact Inspector Table
html = idx_path.read_text(encoding="utf-8")

old_nis = """      let nisAp=0, nisCh=0, nisTot=0;
      for (const r of [...resA, ...resB]) {
        if (intfA.has(r.key) || intfB.has(r.key)) continue;"""
new_nis = """      let nisAp=0, nisCh=0, nisTot=0;
      for (const r of [...resA, ...resB]) {
        if (hasOrdCore && platConf(r) < 0.35) continue;
        if (intfA.has(r.key) || intfB.has(r.key)) continue;"""
assert old_nis in html, "Could not locate NIS loop in index.html"
html = html.replace(old_nis, new_nis, 1)

old_tags = """        if (c.pi > 0) tags.push(`π:${c.pi}`);
        if (!tags.length) tags.push("vdW");"""
new_tags = """        if (c.pi > 0) tags.push(`π:${c.pi}`);
        if (c.ss > 0) tags.push(`SS:${c.ss}`);
        if (c.clash > 0) tags.push(`Clash:${c.clash}`);
        if (!tags.length) tags.push("vdW");"""
assert old_tags in html, "Could not locate contact tags block in index.html"
html = html.replace(old_tags, new_tags, 1)

idx_path.write_text(html, encoding="utf-8")
print("[PASS] Locked IDR-immune NIS & SS/Clash tags in index.html.")

# 2. Patch unippi_cli.py: Remove legacy <50 aa PRODIGY key-mapping guard
cli = cli_path.read_text(encoding="utf-8")
old_guard = """            if min(len(residues_a), len(residues_b)) >= 50:
                if pair_cls == "APOLAR_CHARGED":
                    pair_cls = "CHARGED_APOLAR"
                elif pair_cls == "APOLAR_POLAR":
                    pair_cls = "POLAR_APOLAR\""""
new_guard = """            if pair_cls == "APOLAR_CHARGED":
                pair_cls = "CHARGED_APOLAR"
            elif pair_cls == "APOLAR_POLAR":
                pair_cls = "POLAR_APOLAR\""""
assert old_guard in cli, "Could not locate <50 aa PRODIGY guard in unippi_cli.py"
cli = cli.replace(old_guard, new_guard, 1)
ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] Removed <50 aa PRODIGY guard in unippi_cli.py (ast.parse OK).\n")

# 3. Re-run Verification Tests A & B
zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("="*95)
    print(f"TEST A (FINAL LOCK): unippi_cli.py -i {zip_cands[0].name}")
    print("="*95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)

print("\n" + "="*95)
print("TEST B (FINAL LOCK): unippi_cli.py --generate-demo")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
