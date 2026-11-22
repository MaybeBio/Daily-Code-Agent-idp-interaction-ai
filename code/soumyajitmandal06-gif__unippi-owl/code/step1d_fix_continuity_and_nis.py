#!/usr/bin/env python3
"""
Step 1D: Fixes the 2 mathematical flaws caught by Test B in unippi_cli.py:
  1. Replaces the discontinuous 38% cliff in thermo_gate at -5.5 kcal/mol with
     a smooth continuous ramp: clamp((-dg_uni + 0.5) / 6.0, 0.10, 1.0).
  2. Gates PRODIGY blending by nis_reliability = clamp(nis_total / 20.0, 0.0, 1.0)
     so short peptides/domains without a real non-interacting surface (nis_total < 20)
     rely on 8-channel dG_phys instead of an undefined %NIS fallback, while leaving
     full-length proteins (like fold_vip1_bzip61.zip) 100% unchanged.
"""

import ast
import json
import re
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
src = cli_path.read_text(encoding="utf-8")

# 1. Inspect and fix discontinuous thermo_gate
tg_hits = [(i+1, l) for i, l in enumerate(src.splitlines()) if "thermo_gate" in l]
print("[AUDIT] thermo_gate lines BEFORE fix:", tg_hits)

new_src, n_tg = re.subn(
    r'thermo_gate\s*=\s*clamp\([^\n]+\)\s+if\s+dg_uni\s*>\s*-5\.5\s+else\s+1\.0',
    'thermo_gate = clamp((-dg_uni + 0.5) / 6.0, 0.10, 1.0)',
    src
)
if n_tg == 0:
    # Handle alternate ordering if written as (1.0 if dg_uni <= -5.5 else clamp(...))
    new_src, n_tg = re.subn(
        r'thermo_gate\s*=\s*1\.0\s+if\s+dg_uni\s*<=\s*-5\.5\s+else\s+clamp\([^\n]+\)',
        'thermo_gate = clamp((-dg_uni + 0.5) / 6.0, 0.10, 1.0)',
        src
    )
assert n_tg >= 1, "Could not locate discontinuous thermo_gate line in unippi_cli.py!"
print(f"[OK] Replaced discontinuous thermo_gate with smooth continuous ramp (matches={n_tg}).")

# 2. Inspect and add nis_reliability gating to PRODIGY blend
gate_hits = [(i+1, l) for i, l in enumerate(new_src.splitlines()) if "prodigy_gate" in l or "nis_tot" in l or "nis_total" in l]
print("[AUDIT] NIS & gate lines BEFORE fix:", gate_hits[:8])

# Detect exact variable name for NIS count (nis_total vs nis_tot)
nis_var = "nis_total" if "nis_total" in new_src else ("nis_tot" if "nis_tot" in new_src else None)
assert nis_var is not None, "Could not find nis_total / nis_tot variable in unippi_cli.py!"

# Multiply prodigy_gate (or contact_gate) by nis_reliability = clamp(nis_var / 20.0, 0.0, 1.0)
gate_var = None
for cand in ("prodigy_gate", "contact_gate", "ic_gate"):
    if f"{cand} =" in new_src:
        gate_var = cand
        break
assert gate_var is not None, "Could not find prodigy_gate / contact_gate variable in unippi_cli.py!"

if "nis_reliability" not in new_src:
    lines = new_src.splitlines()
    idx_g = next(i for i, l in enumerate(lines) if f"{gate_var} =" in l)
    ind = lines[idx_g][:len(lines[idx_g]) - len(lines[idx_g].lstrip())]
    # Insert nis_reliability right after the gate assignment block
    insert_at = idx_g + 1
    while insert_at < len(lines) and (lines[insert_at].strip().startswith("else:") or lines[insert_at].strip().startswith(f"{gate_var} =")):
        insert_at += 1
    lines.insert(insert_at, f"{ind}nis_reliability = clamp({nis_var} / 20.0, 0.0, 1.0)")
    lines.insert(insert_at + 1, f"{ind}{gate_var} = {gate_var} * nis_reliability")
    new_src = "\n".join(lines) + "\n"
    print(f"[OK] Added nis_reliability = clamp({nis_var} / 20.0, 0.0, 1.0) to {gate_var}.")

# 3. Validate syntax with ast.parse() before writing
ast.parse(new_src)
cli_path.write_text(new_src, encoding="utf-8")
print("[PASS] unippi_cli.py updated and verified cleanly via ast.parse().\n")

# 4. Run Test A & Test B and print full per-model diagnostic tables
def print_audit_table(json_path):
    if not json_path.exists():
        return
    data = json.loads(json_path.read_text(encoding="utf-8"))
    print(f"\n[Audit Table: {json_path.parent.name}]")
    print(f"{'Model':<38} | {'PRODIGY':>7} | {'dG_phys':>7} | {'dG_uni':>7} | {'Kd (nM)':>10} | {'wBSA':>7} | {'Pro':>3} | {'ipTM':>5} | {'LIS':>5} | {'UIS':>5} | {'Pattern':<9}")
    print("-"*128)
    for r in data:
        kd_s = f"{r['Kd_nM']:.2f}" if r['Kd_nM'] < 1e7 else ">10mM(NB)"
        print(f"{r['complex_id'][:38]:<38} | {r['dG_prodigy_kcal']:>7.2f} | {r['dG_phys_kcal']:>7.2f} | {r['dG_uni_kcal']:>7.2f} | {kd_s:>10} | {r['wBSA_A2']:>7.1f} | {r['helical_prolines']:>3d} | {r['ipTM']:>5.3f} | {r['LIS']:>5.3f} | {r['UIS']:>5.3f} | {r['pattern_id']:<9}")

zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("="*95)
    print(f"TEST A: REAL 5-MODEL AF3 ARCHIVE ({zip_cands[0].name})")
    print("="*95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)
    print_audit_table(Path("./af3_vip1_bzip61_verified/unippi_full_report.json"))

print("\n" + "="*95)
print("TEST B: POSITIVE vs. NEGATIVE CONTROL BENCHMARK (--generate-demo)")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
print_audit_table(Path("./results_demo_verified/unippi_full_report.json"))
