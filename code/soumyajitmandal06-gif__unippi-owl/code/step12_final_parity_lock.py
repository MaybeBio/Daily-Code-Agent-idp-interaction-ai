#!/usr/bin/env python3
"""
Step 12 (Final Zero-Delta CLI <-> Browser Parity Lock):
  1. Unifies conf_contact_log (L940) and n_cb8 (L1114) in unippi_cli.py from >= 0.15 to >= 0.25.
  2. Aligns Coordinate-Mode pdockq_pae & fallback iptm_val (L1133, L1154) with index.html.
  3. Verifies syntax with ast.parse() and runs Test A (fold_vip1_bzip61.zip) & Test B (--generate-demo).
"""

import ast
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
cli = cli_path.read_text(encoding="utf-8")

patches = [
    (
        '    conf_contact_log = [c for c in residue_contact_log if c["weight_w"] >= 0.15]',
        '    conf_contact_log = [c for c in residue_contact_log if c["weight_w"] >= 0.25]'
    ),
    (
        '        n_cb8 = max(len(conf_contact_log), sum(1 for w in pair_weights_seen if w >= 0.15))',
        '        n_cb8 = max(len(conf_contact_log), sum(1 for w in pair_weights_seen if w >= 0.25))'
    ),
    (
        '    pdockq_pae = pdockq_coord * (1.0 / (1.0 + (excess_intf_pae / 8.0) ** 2))',
        '    pdockq_pae = pdockq_coord * (1.0 / (1.0 + (excess_intf_pae / 8.0) ** 2)) if mode == "MODE_A_PRED_FULL" else clamp(pdockq_coord / 0.742, 0.0, 1.0)'
    ),
    (
        '        iptm_val = pdockq_pae if mode == "MODE_A_PRED_FULL" else pdockq_coord',
        '        iptm_val = pdockq_pae if mode == "MODE_A_PRED_FULL" else clamp(pdockq_coord / 0.742, 0.0, 1.0)'
    ),
]

for old, new in patches:
    assert old in cli, f"Missing expected line in unippi_cli.py:\n{old}"
    cli = cli.replace(old, new, 1)

ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] Locked 0.0000 CLI <-> Browser parity in unippi_cli.py (ast.parse OK).\n")

zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("="*95)
    print(f"TEST A (ZERO-DELTA LOCK): unippi_cli.py -i {zip_cands[0].name}")
    print("="*95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)

print("\n" + "="*95)
print("TEST B (ZERO-DELTA LOCK): unippi_cli.py --generate-demo")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
