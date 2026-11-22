#!/usr/bin/env python3
"""
AST-introspective repair for NameError ('is_residue_helical') in unippi_cli.py:
  1. Discovers the actual helix helper function (name + arity) from unippi_cli.py's AST.
  2. Patches line 927 to call the verified symbol (or inserts a self-contained math.dist helper).
  3. Prints the actual thermodynamics & blending code around 15.5 * math.tanh.
  4. Runs Test A (fold_vip1_bzip61.zip) and Test B (--generate-demo).
"""

import ast
import json
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
src = cli_path.read_text(encoding="utf-8")
tree = ast.parse(src)

# 1. Introspect all FunctionDef nodes in unippi_cli.py for helix detection
hel_funcs = [
    (node.name, [a.arg for a in node.args.args], node.lineno)
    for node in ast.walk(tree)
    if isinstance(node, ast.FunctionDef) and "hel" in node.name.lower()
]
print("[AST AUDIT] Helix-related functions found in unippi_cli.py:", hel_funcs)

if hel_funcs:
    real_name, real_args, real_line = hel_funcs[0]
    if len(real_args) == 2:
        new_call = f"{real_name}(idx, chain_list)"
    else:
        new_call = f"{real_name}(idx, chain_list, mode)"
    src = src.replace("is_residue_helical(idx, chain_list, mode)", new_call)
    print(f"[OK] Replaced undefined 'is_residue_helical(...)' with verified '{new_call}' (defined at L{real_line}).")
else:
    # Self-contained fallback using math.dist on Residue centroids (checks 3.6-res alpha-helix turn)
    lines = src.splitlines()
    idx_pro = next(i for i, l in enumerate(lines) if "is_residue_helical(" in l)
    ind = lines[idx_pro][:len(lines[idx_pro]) - len(lines[idx_pro].lstrip())]
    # Find start of helical_prolines = 0 to define a clean local helper right above the loop
    idx_start = idx_pro
    while idx_start > 0 and "helical_prolines = 0" not in lines[idx_start]:
        idx_start -= 1
    b_ind = lines[idx_start][:len(lines[idx_start]) - len(lines[idx_start].lstrip())]
    helper_lines = [
        f'{b_ind}def _is_interior_helix(i_res, c_list):',
        f'{b_ind}    if i_res < 4 or i_res >= len(c_list) - 4:',
        f'{b_ind}        return False',
        f'{b_ind}    c0 = getattr(c_list[i_res], "centroid", None) or getattr(c_list[i_res], "ca_coord", None)',
        f'{b_ind}    cp = getattr(c_list[i_res - 4], "centroid", None) or getattr(c_list[i_res - 4], "ca_coord", None)',
        f'{b_ind}    cn = getattr(c_list[i_res + 4], "centroid", None) or getattr(c_list[i_res + 4], "ca_coord", None)',
        f'{b_ind}    if c0 is None or cp is None or cn is None:',
        f'{b_ind}        return False',
        f'{b_ind}    d_prev, d_next, d_span = math.dist(cp, c0), math.dist(c0, cn), math.dist(cp, cn)',
        f'{b_ind}    return (5.0 <= d_prev <= 7.0) and (5.0 <= d_next <= 7.0) and (d_span >= 10.2)',
    ]
    lines[idx_pro] = lines[idx_pro].replace("is_residue_helical(idx, chain_list, mode)", "_is_interior_helix(idx, chain_list)")
    lines = lines[:idx_start] + helper_lines + lines[idx_start:]
    src = "\n".join(lines) + "\n"
    print("[OK] Injected self-contained '_is_interior_helix(idx, chain_list)' using math.dist.")

# Verify syntax AND zero remaining references to is_residue_helical
assert "is_residue_helical" not in src, "FAIL: is_residue_helical still present!"
ast.parse(src)
cli_path.write_text(src, encoding="utf-8")
print("[PASS] Symbol fix verified via AST.\n")

# 2. Print the exact thermodynamics & blending block around 15.5 * math.tanh
lines = src.splitlines()
idx_tanh = next((i for i, l in enumerate(lines) if "math.tanh" in l), -1)
if idx_tanh >= 0:
    print("="*95)
    print(f"CODE AUDIT: THERMODYNAMICS & BLENDING (Lines {max(1, idx_tanh-15)}..{min(len(lines), idx_tanh+30)})")
    print("="*95)
    for k in range(max(0, idx_tanh - 15), min(len(lines), idx_tanh + 30)):
        print(f"  L{k+1:4d}: {lines[k]}")

# 3. Run Test A (Real 5-model AF3 zip) & Test B (Positive/Negative Controls)
zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("\n" + "="*95)
    print(f"TEST A: {zip_cands[0].name}")
    print("="*95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)
    rep_path = Path("./af3_vip1_bzip61_verified/unippi_full_report.json")
    if rep_path.exists():
        data = json.loads(rep_path.read_text(encoding="utf-8"))
        print("\n[Per-Model Audit Table across all 5 AF3 Models]")
        print(f"{'Model':<26} | {'dG_PRODIGY':>10} | {'dG_phys':>8} | {'dG_uni':>7} | {'wBSA':>7} | {'SB':>3} | {'HB':>3} | {'HP':>3} | {'Pro':>3} | {'ipTM':>5} | {'LIS':>5} | {'UIS':>5} | {'Pattern':<9}")
        print("-"*122)
        for r in data:
            print(f"{r['complex_id']:<26} | {r['dG_prodigy_kcal']:>10.2f} | {r['dG_phys_kcal']:>8.2f} | {r['dG_uni_kcal']:>7.2f} | {r['wBSA_A2']:>7.1f} | {r['salt_bridges']:>3d} | {r['h_bonds']:>3d} | {r['hydrophobic_contacts']:>3d} | {r['helical_prolines']:>3d} | {r['ipTM']:>5.3f} | {r['LIS']:>5.3f} | {r['UIS']:>5.3f} | {r['pattern_id']:<9}")

print("\n" + "="*95)
print("TEST B: POSITIVE vs. NEGATIVE CONTROLS (--generate-demo)")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
