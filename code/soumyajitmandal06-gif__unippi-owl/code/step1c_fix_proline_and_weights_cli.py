#!/usr/bin/env python3
"""
Step 1C:
  1. Audits and fixes helical_prolines (strict interior helix + confident interface +
     strain-gated penalty so already-relaxed AF3 Prolines at pLDDT>=85 & PAE<=3.5 A
     are never double-penalized).
  2. Locks phys_weight cleanly in [0.50, 0.70] with zero stray overrides.
  3. Re-runs Test A (fold_vip1_bzip61.zip) and Test B (--generate-demo).
"""

import ast
import json
import re
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
cli = cli_path.read_text(encoding="utf-8")

# 1. Inspect existing helical_prolines & phys_weight lines for transparency
for tag in ("helical_prolines", "phys_weight"):
    hits = [f"L{i+1}: {l}" for i, l in enumerate(cli.splitlines()) if tag in l]
    print(f"[AUDIT] {tag} occurrences ({len(hits)}):")
    for h in hits[:8]:
        print("  ", h)

# 2. Tighten is_residue_helical so helix caps / turns aren't falsely flagged as interior helix kinks
old_hel = """    if idx + 4 < len(chain_residues) and idx - 4 >= 0:
        d_next = vec_dist(chain_residues[idx + 4].centroid, ci)
        d_prev = vec_dist(chain_residues[idx - 4].centroid, ci)
        if 4.8 <= d_next <= 7.2 and 4.8 <= d_prev <= 7.2:
            return True"""
new_hel = """    if idx + 4 < len(chain_residues) and idx - 4 >= 0:
        d_next = vec_dist(chain_residues[idx + 4].centroid, ci)
        d_prev = vec_dist(chain_residues[idx - 4].centroid, ci)
        d_span = vec_dist(chain_residues[idx - 4].centroid, chain_residues[idx + 4].centroid)
        if 5.0 <= d_next <= 7.0 and 5.0 <= d_prev <= 7.0 and d_span >= 10.2:
            return True"""
if old_hel in cli:
    cli = cli.replace(old_hel, new_hel)

# 3. Make helical_prolines count ONLY strained/unrelaxed interior Prolines near the interface
#    (pLDDT < 85 or mean pair PAE > 3.5 A), avoiding double-penalizing relaxed AF3 models
old_pro_loop = re.compile(
    r'    helical_prolines = 0\n.*?helical_prolines \+= 1',
    re.DOTALL
)
new_pro_loop = """    helical_prolines = 0
    _mean_pae_tmp = mean_list(intf_pae_vals, 3.5)
    for chain_list, intf_set in ((residues_a, intf_res_a), (residues_b, intf_res_b)):
        for idx, res in enumerate(chain_list):
            if res.canonical_name == "PRO" and residue_plateau_conf(res) >= 0.50 and is_residue_helical(idx, chain_list, mode):
                near = any(abs(res.res_seq - int(k.split("|")[1])) <= 2 for k in intf_set)
                is_strained = (res.mean_bfactor < 85.0) or (_mean_pae_tmp > 3.5)
                if near and is_strained:
                    helical_prolines += 1"""
cli, n_pro = old_pro_loop.subn(new_pro_loop, cli, count=1)
print(f"\n[OK] Updated helical_prolines strain-gating block (matches={n_pro})")

# 4. Lock phys_weight cleanly to clamp(0.50 + 0.20 * apolar_frac, 0.50, 0.70) (no stray 0.82 override)
cli = re.sub(
    r'phys_weight\s*=\s*clamp\([^\n]+\)',
    'phys_weight = clamp(0.50 + 0.20 * apolar_frac + (0.10 if helical_prolines > 0 else 0.0), 0.50, 0.70)',
    cli
)
# Also remove any separate phys_weight = min(0.95, ...) or w_bsa bump
cli = re.sub(
    r'if\s+w_bsa\s*>=\s*800\.0:\s*\n\s*phys_weight\s*=[^\n]+',
    'if False:\n        pass',
    cli
)

# 5. Validate with ast.parse before saving
ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] unippi_cli.py syntax verified via ast.parse().\n")

# 6. Run Test A & Test B
zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("="*95)
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
