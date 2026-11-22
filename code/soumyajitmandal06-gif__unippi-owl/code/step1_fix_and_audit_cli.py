#!/usr/bin/env python3
"""
Step 1 of 2: Fix Ordered-Core SASA Parity, Unified w_min = 0.25, and
Kd formatting in unippi_cli.py, then run a full numerical audit on
Prof. Daisuke's real 5-model AF3 archive (fold_vip1_bzip61.zip) + demo controls.
"""

import re
import subprocess
import sys
from pathlib import Path

cli_path = Path("unippi_cli.py")
cli = cli_path.read_text(encoding="utf-8")

# 1. Flexibly patch Ordered-Core SASA & unified w_min = 0.25 in analyze_interface_pair
sasa_pattern = re.compile(
    r'    total_free_a = sum\([^\n]+\)\n'
    r'    total_free_b = sum\([^\n]+\)\n'
    r'    total_comp_ab = sum\([^\n]+\)\n'
    r'    bsa = max\(0\.0,[^\n]+\)\n'
    r'    mean_w_intf = [^\n]+\n'
    r'    w_bsa = [^\n]+',
    re.MULTILINE
)

ordered_sasa_block = """    conf_contacts_core = [c for c in residue_contact_log if c["weight_w"] >= 0.25]
    has_ord_core = len(conf_contacts_core) > 0
    ord_keys_a = {r.key for r in residues_a if residue_plateau_conf(r) >= 0.35} if has_ord_core else {r.key for r in residues_a}
    ord_keys_b = {r.key for r in residues_b if residue_plateau_conf(r) >= 0.35} if has_ord_core else {r.key for r in residues_b}
    total_free_a = sum(v for k, v in sasa_free_a.items() if k in ord_keys_a)
    total_free_b = sum(v for k, v in sasa_free_b.items() if k in ord_keys_b)
    total_comp_ab = sum(sasa_complex.get(r.key, 0.0) for r in (residues_a + residues_b) if (r.key in ord_keys_a or r.key in ord_keys_b))
    bsa = max(0.0, (total_free_a + total_free_b) - total_comp_ab)
    mean_w_intf = mean_list([c["weight_w"] for c in conf_contacts_core], 0.0) if has_ord_core else mean_list(pair_weights_seen, 0.0)
    w_bsa = mean_w_intf * bsa"""

if "ord_keys_a =" not in cli:
    new_cli, n_sub = sasa_pattern.subn(ordered_sasa_block, cli, count=1)
    if n_sub == 0:
        # Fallback: find lines around 'bsa = max(0.0'
        lines = cli.splitlines()
        idx_bsa = next((i for i, l in enumerate(lines) if "bsa = max(0.0," in l), -1)
        if idx_bsa < 0:
            raise RuntimeError("Could not locate bsa = max(0.0, ...) in unippi_cli.py")
        # Print surrounding context and replace block
        start_i = max(0, idx_bsa - 4)
        end_i = min(len(lines), idx_bsa + 3)
        print("[DEBUG] Found BSA block:\n" + "\n".join(lines[start_i:end_i]))
        while start_i < idx_bsa and "free_a" not in lines[start_i]:
            start_i += 1
        while end_i > idx_bsa and "w_bsa" not in lines[end_i - 1]:
            end_i -= 1
        lines[start_i:end_i] = ordered_sasa_block.splitlines()
        cli = "\n".join(lines) + "\n"
    else:
        cli = new_cli
    print("[OK] Patched ordered-core SASA block in unippi_cli.py")
else:
    cli = cli.replace('c["weight_w"] >= 0.35', 'c["weight_w"] >= 0.25')
    print("[OK] Updated existing ordered-core SASA block to unified w_min = 0.25")

# 2. Unify w_min = 0.25 for intf_residues and Step 1 locus printout + compute real pct_disordered
cli = cli.replace('c["weight_w"] >= 0.35', 'c["weight_w"] >= 0.25')
cli = cli.replace('c.get("weight_w", 0) >= 0.35', 'c.get("weight_w", 0) >= 0.25')
cli = cli.replace('Ordered Core Loci (w>=0.35):', 'Ordered Core Loci (w>=0.25):')

# Add pct_disordered calculation right next to len_a/len_b in the returned record
if '"pct_disordered":' not in cli:
    cli = cli.replace(
        '"len_a": len(residues_a), "len_b": len(residues_b),',
        '"len_a": len(residues_a), "len_b": len(residues_b), "pct_disordered": round(100.0 * sum(1 for r in (residues_a + residues_b) if r.mean_bfactor < 70.0) / max(1, len(residues_a) + len(residues_b)), 1),'
    )

# 3. Fix Kd formatting in CLI table & Step 4 printout (show pM when < 0.05 nM instead of 0.00 nM)
cli = cli.replace(
    'f"{r[\'Kd_nM\']:<11.2f}" if r[\'Kd_nM\'] < 1e7 else ">10mM (NB) "',
    'f"{r[\'Kd_nM\']*1000:.1f} pM".ljust(11) if r[\'Kd_nM\'] < 0.05 else (f"{r[\'Kd_nM\']:<11.2f}" if r[\'Kd_nM\'] < 1e7 else ">10mM (NB) ")'
)
cli = cli.replace(
    "(Kd = {top['Kd_nM']:.2f} nM)",
    "(Kd = {f\"{top['Kd_nM']*1000:.2f} pM\" if top['Kd_nM'] < 0.05 else f\"{top['Kd_nM']:.2f} nM\"})"
)

cli_path.write_text(cli, encoding="utf-8")

# Verify critical markers in unippi_cli.py
c_now = cli_path.read_text(encoding="utf-8")
assert "ord_keys_a =" in c_now, "FAIL: ord_keys_a missing"
assert "16.5 * math.tanh(e_raw / 16.5)" in c_now, "FAIL: tanh compensation missing"
assert 'c["weight_w"] >= 0.25' in c_now, "FAIL: unified w_min = 0.25 missing"
print("[PASS] unippi_cli.py passed all parity & physics assertions.\n")

# 4. Run Live Verification #1: Prof. Daisuke's real 5-model AF3 zip archive
zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print("="*95)
    print(f"TEST A: REAL 5-MODEL AF3 ARCHIVE ({zip_cands[0].name})")
    print("="*95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)
else:
    print("[WARN] fold_vip1_bzip61.zip not found in current directory.")

# 5. Run Live Verification #2: Positive & Negative Control Benchmark
print("\n" + "="*95)
print("TEST B: POSITIVE vs. NEGATIVE CONTROL BENCHMARK (--generate-demo)")
print("="*95)
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
