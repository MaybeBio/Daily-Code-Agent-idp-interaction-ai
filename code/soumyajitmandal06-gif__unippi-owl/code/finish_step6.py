#!/usr/bin/env python3
"""
Completes Step 6 cleanly (idempotent):
  1. Verifies Part 1 patches are already present in index.html.
  2. Skips blank lines above 'w_sb +=' in unippi_cli.py and gates the integer
     contact counters (n_sb..n_clash) under 'if w_ij >= 0.25:'.
  3. Validates via ast.parse() and runs Test A (fold_vip1_bzip61.zip) & Test B (--generate-demo).
"""

import ast
import re
import subprocess
import sys
from pathlib import Path

idx_path = Path("index.html")
cli_path = Path("unippi_cli.py")

# 1. Verify index.html already has all 3 Step 6 upgrades
h = idx_path.read_text(encoding="utf-8")
assert "if (wij >= 0.25) { nSb += effSb;" in h, "index.html missing wij >= 0.25 count gate"
assert 'id="hideLowConfContactsCheck"' in h, "index.html missing Ordered-only contact toggle"
assert "wHP=${(r.wHp||0).toFixed(1)}" in h, "index.html missing wHP in Step 4"
print("[PASS] Verified all 3 Step 6 upgrades are active in index.html.")

# 2. Patch unippi_cli.py (skip blank lines above w_sb +=)
cli_lines = cli_path.read_text(encoding="utf-8").splitlines()
idx_wsb = next((i for i, l in enumerate(cli_lines) if "w_sb +=" in l), -1)
assert idx_wsb > 0, "Could not locate 'w_sb +=' in unippi_cli.py"

m_wvar = re.search(r'w_sb\s*\+=\s*([A-Za-z_][A-Za-z0-9_]*)\s*\*', cli_lines[idx_wsb])
assert m_wvar, f"Could not parse weight variable from: {cli_lines[idx_wsb]}"
w_var = m_wvar.group(1)

# Walk past any blank lines directly above w_sb +=
cnt_end = idx_wsb
while cnt_end > 0 and not cli_lines[cnt_end - 1].strip():
    cnt_end -= 1

# Walk backward across all contiguous integer counter lines (n_sb..n_clash)
cnt_start = cnt_end
while cnt_start > 0 and re.match(r'^\s*n_(?:sb|hb|hp|pi|ss|clash|salt|hbond|hydro)\s*\+=', cli_lines[cnt_start - 1]):
    cnt_start -= 1

assert cnt_start < cnt_end, f"No integer counter lines found above L{idx_wsb+1}: {cli_lines[max(0,cnt_end-6):idx_wsb+1]}"

if f"if {w_var} >= 0.25:" not in cli_lines[cnt_start - 1]:
    ind = cli_lines[cnt_start][:len(cli_lines[cnt_start]) - len(cli_lines[cnt_start].lstrip())]
    print(f"[INFO] Gating {cnt_end - cnt_start} integer counter lines (L{cnt_start+1}..L{cnt_end}) under 'if {w_var} >= 0.25:':")
    for l in cli_lines[cnt_start:cnt_end]:
        print("  ", l)
    gated = [f"{ind}if {w_var} >= 0.25:"] + [f"    {l}" for l in cli_lines[cnt_start:cnt_end]]
    cli_lines = cli_lines[:cnt_start] + gated + cli_lines[cnt_end:]
    cli_new = "\n".join(cli_lines) + "\n"
    ast.parse(cli_new)
    cli_path.write_text(cli_new, encoding="utf-8")
    print(f"[PASS] Patched unippi_cli.py cleanly (ast.parse OK).\n")
else:
    ast.parse("\n".join(cli_lines) + "\n")
    print("[INFO] Integer counters already gated in unippi_cli.py (ast.parse OK).\n")

# 3. Run Live Verification Tests A & B
zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("="*95)
    print(f"TEST A (LOCKED BASELINE): unippi_cli.py -i {zip_cands[0].name}")
    print("="*95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)

print("\n" + "="*95)
print("TEST B (LOCKED BASELINE): unippi_cli.py --generate-demo")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
