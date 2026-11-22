#!/usr/bin/env python3
"""
Indentation-safe repair of the helical_prolines block in unippi_cli.py
+ full printout of the dg_prodigy / phys_weight / dg_uni equations
+ live verification on fold_vip1_bzip61.zip (Test A) and --generate-demo (Test B).
"""

import ast
import json
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
lines = cli_path.read_text(encoding="utf-8").splitlines()

# 1. Locate the exact helical_prolines block by scanning for intf_pae_vals or helical_prolines = 0
idx_err = next((i for i, l in enumerate(lines) if "intf_pae_vals" in l), -1)
if idx_err >= 0:
    start_i = idx_err
    while start_i > 0 and "helical_prolines = 0" not in lines[start_i]:
        start_i -= 1
    end_i = idx_err
    while end_i < len(lines) and "helical_prolines += 1" not in lines[end_i]:
        end_i += 1
    end_i += 1  # inclusive of helical_prolines += 1

    print(f"[INFO] Found broken helical_prolines block at lines {start_i+1}..{end_i}:")
    for k in range(start_i, end_i):
        print(f"  L{k+1:4d}: {lines[k]}")

    ind = lines[start_i][:len(lines[start_i]) - len(lines[start_i].lstrip())]
    replacement = [
        f'{ind}helical_prolines = 0',
        f'{ind}_conf_paes_pro = [c.get("pae_A", c.get("pae", c.get("mean_pae", 3.5))) for c in residue_contact_log if c.get("weight_w", 0.0) >= 0.25]',
        f'{ind}_mean_pae_tmp = mean_list(_conf_paes_pro, mean_list(intf_pae_values, 3.5))',
        f'{ind}_ord_intf_a = {{r.res_seq for r in residues_a if r.key in intf_res_a and residue_plateau_conf(r) >= 0.50}}',
        f'{ind}_ord_intf_b = {{r.res_seq for r in residues_b if r.key in intf_res_b and residue_plateau_conf(r) >= 0.50}}',
        f'{ind}for chain_list, ord_seqs in ((residues_a, _ord_intf_a), (residues_b, _ord_intf_b)):',
        f'{ind}    for idx, res in enumerate(chain_list):',
        f'{ind}        if res.canonical_name == "PRO" and residue_plateau_conf(res) >= 0.50 and is_residue_helical(idx, chain_list, mode):',
        f'{ind}            near = any(abs(res.res_seq - s_num) <= 2 for s_num in ord_seqs)',
        f'{ind}            is_strained = (res.mean_bfactor < 80.0) or (_mean_pae_tmp > 4.0)',
        f'{ind}            if near and is_strained:',
        f'{ind}                helical_prolines += 1',
    ]
    lines = lines[:start_i] + replacement + lines[end_i:]
    new_src = "\n".join(lines) + "\n"
    ast.parse(new_src)
    cli_path.write_text(new_src, encoding="utf-8")
    print(f"[PASS] Replaced lines {start_i+1}..{end_i} (indent={len(ind)} spaces) & verified via ast.parse().\n")
else:
    print("[INFO] 'intf_pae_vals' not found; validating current unippi_cli.py via ast.parse()...")
    ast.parse("\n".join(lines) + "\n")

# 2. Print the exact thermodynamics & blending block so we solve the 81.7% weight mystery with zero guesswork
idx_prod = next((i for i, l in enumerate(lines) if "dg_prodigy =" in l), -1)
if idx_prod >= 0:
    print("="*95)
    print(f"CODE AUDIT: THERMODYNAMICS & BLENDING BLOCK (Lines {idx_prod+1}..{min(len(lines), idx_prod+45)})")
    print("="*95)
    for k in range(idx_prod, min(len(lines), idx_prod + 45)):
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
