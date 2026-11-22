#!/usr/bin/env python3
"""
Repairs the SASA block indentation in unippi_cli.py with ast.parse() verification,
then runs Test A (Prof. Daisuke's fold_vip1_bzip61.zip) and Test B (--generate-demo).
"""

import ast
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
lines = cli_path.read_text(encoding="utf-8").splitlines()

# 1. Inspect lines around the inserted conf_contacts_core / line 939
idx_ins = next((i for i, l in enumerate(lines) if "conf_contacts_core =" in l), -1)
if idx_ins < 0:
    print("[ERROR] Could not find conf_contacts_core in unippi_cli.py. Showing lines 920-955:")
    for n in range(max(0, 919), min(len(lines), 955)):
        print(f"{n+1:4d}: {lines[n]}")
    sys.exit(1)

print("[INFO] Broken region BEFORE fix:")
for n in range(max(0, idx_ins - 3), min(len(lines), idx_ins + 18)):
    print(f"{n+1:4d}: {lines[n]}")

# 2. Detect true surrounding indentation from the original 'bsa = max(0.0,' line below our insert
idx_orig_bsa = next((i for i in range(idx_ins + 8, min(len(lines), idx_ins + 20)) if "bsa = max(0.0," in lines[i]), -1)
if idx_orig_bsa >= 0:
    true_indent = lines[idx_orig_bsa][:len(lines[idx_orig_bsa]) - len(lines[idx_orig_bsa].lstrip())]
    # Find end of the original BSA/w_bsa block (up to 6 lines after original bsa)
    end_replace = idx_orig_bsa + 1
    for k in range(idx_orig_bsa + 1, min(len(lines), idx_orig_bsa + 7)):
        if any(tok in lines[k] for tok in ("mean_w_intf", "w_bsa", "wbsa")):
            end_replace = k + 1
else:
    # Fallback: take indentation from previous non-empty, non-comment code line
    prev_code = next((lines[i] for i in range(idx_ins - 1, -1, -1) if lines[i].strip() and not lines[i].lstrip().startswith("#")), "    ")
    true_indent = prev_code[:len(prev_code) - len(prev_code.lstrip())]
    end_replace = idx_ins + 9

# Also check if any old 'total_free_a' / 'total_comp_ab' lines remain right before idx_ins
start_replace = idx_ins
while start_replace > 0 and any(tok in lines[start_replace - 1] for tok in ("total_free_a", "total_free_b", "total_comp_ab")):
    start_replace -= 1

clean_sasa_lines = [
    f'{true_indent}conf_contacts_core = [c for c in residue_contact_log if c["weight_w"] >= 0.25]',
    f'{true_indent}has_ord_core = len(conf_contacts_core) > 0',
    f'{true_indent}ord_keys_a = {{r.key for r in residues_a if residue_plateau_conf(r) >= 0.35}} if has_ord_core else {{r.key for r in residues_a}}',
    f'{true_indent}ord_keys_b = {{r.key for r in residues_b if residue_plateau_conf(r) >= 0.35}} if has_ord_core else {{r.key for r in residues_b}}',
    f'{true_indent}total_free_a = sum(v for k, v in sasa_free_a.items() if k in ord_keys_a)',
    f'{true_indent}total_free_b = sum(v for k, v in sasa_free_b.items() if k in ord_keys_b)',
    f'{true_indent}total_comp_ab = sum(sasa_complex.get(r.key, 0.0) for r in (residues_a + residues_b) if (r.key in ord_keys_a or r.key in ord_keys_b))',
    f'{true_indent}bsa = max(0.0, (total_free_a + total_free_b) - total_comp_ab)',
    f'{true_indent}mean_w_intf = mean_list([c["weight_w"] for c in conf_contacts_core], 0.0) if has_ord_core else mean_list(pair_weights_seen, 0.0)',
    f'{true_indent}w_bsa = mean_w_intf * bsa',
]

new_lines = lines[:start_replace] + clean_sasa_lines + lines[end_replace:]
new_source = "\n".join(new_lines) + "\n"

# 3. Validate syntax with ast.parse BEFORE writing
try:
    ast.parse(new_source)
except SyntaxError as e:
    print(f"\n[FAIL] SyntaxError after repair at line {e.lineno}: {e.msg}")
    for n in range(max(0, (e.lineno or 1) - 6), min(len(new_lines), (e.lineno or 1) + 6)):
        print(f"{n+1:4d}: {new_lines[n]}")
    sys.exit(1)

cli_path.write_text(new_source, encoding="utf-8")
print(f"\n[PASS] Repaired unippi_cli.py with indent={len(true_indent)} spaces and verified cleanly via ast.parse()!\n")

# 4. Run Test A (Real 5-model AF3 zip) & Test B (Positive/Negative Controls)
zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("="*95)
    print(f"TEST A: REAL 5-MODEL AF3 ARCHIVE ({zip_cands[0].name})")
    print("="*95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)

print("\n" + "="*95)
print("TEST B: POSITIVE vs. NEGATIVE CONTROL BENCHMARK (--generate-demo)")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
