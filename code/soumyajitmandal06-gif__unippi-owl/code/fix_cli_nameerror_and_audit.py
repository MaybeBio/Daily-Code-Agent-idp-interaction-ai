#!/usr/bin/env python3
"""
Fixes intf_pae_vals -> confident-core PAE (fallback intf_pae_values) in unippi_cli.py,
prints the exact phys_weight/dg_uni block for transparency, and runs Test A & Test B.
"""

import ast
import json
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
cli = cli_path.read_text(encoding="utf-8")

old_block = """    helical_prolines = 0
    _mean_pae_tmp = mean_list(intf_pae_vals, 3.5)
    for chain_list, intf_set in ((residues_a, intf_res_a), (residues_b, intf_res_b)):
        for idx, res in enumerate(chain_list):
            if res.canonical_name == "PRO" and residue_plateau_conf(res) >= 0.50 and is_residue_helical(idx, chain_list, mode):
                near = any(abs(res.res_seq - int(k.split("|")[1])) <= 2 for k in intf_set)
                is_strained = (res.mean_bfactor < 85.0) or (_mean_pae_tmp > 3.5)
                if near and is_strained:
                    helical_prolines += 1"""

new_block = """    helical_prolines = 0
    _conf_paes_pro = [
        c.get("pae_A", c.get("pae", c.get("mean_pae", 3.5)))
        for c in residue_contact_log
        if c.get("weight_w", 0.0) >= 0.25
    ]
    _mean_pae_tmp = mean_list(_conf_paes_pro, mean_list(intf_pae_values, 3.5))
    _ord_intf_a = {r.res_seq for r in residues_a if r.key in intf_res_a and residue_plateau_conf(r) >= 0.50}
    _ord_intf_b = {r.res_seq for r in residues_b if r.key in intf_res_b and residue_plateau_conf(r) >= 0.50}
    for chain_list, ord_seqs in ((residues_a, _ord_intf_a), (residues_b, _ord_intf_b)):
        for idx, res in enumerate(chain_list):
            if res.canonical_name == "PRO" and residue_plateau_conf(res) >= 0.50 and is_residue_helical(idx, chain_list, mode):
                near = any(abs(res.res_seq - s_num) <= 2 for s_num in ord_seqs)
                is_strained = (res.mean_bfactor < 80.0) or (_mean_pae_tmp > 4.0)
                if near and is_strained:
                    helical_prolines += 1"""

assert old_block in cli, "Could not locate exact helical_prolines block in unippi_cli.py"
cli = cli.replace(old_block, new_block, 1)

ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] Fixed NameError & upgraded strain-gating to confident-core PAE (ast.parse OK).\n")

# Audit exact dg_uni / phys_weight lines so we verify zero stray weight bumps
print("[AUDIT] Active phys_weight & dg_uni lines in unippi_cli.py:")
for i, line in enumerate(cli.splitlines(), 1):
    if any(tok in line for tok in ("phys_weight", "dg_uni =", "dg_blended")):
        print(f"  L{i:4d}: {line}")

zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("\n" + "="*95)
    print(f"TEST A (FINAL CLI CHECK): {zip_cands[0].name}")
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
print("TEST B (FINAL CLI CHECK): POSITIVE vs. NEGATIVE CONTROLS (--generate-demo)")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
