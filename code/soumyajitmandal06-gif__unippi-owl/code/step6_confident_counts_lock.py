#!/usr/bin/env python3
"""
Step 6 (Final Phase 1 Hygiene Lock):
  1. Gates displayed integer contact counts (SB, HB, HP, pi, SS, clash) by wij >= 0.25
     in BOTH index.html and unippi_cli.py so disordered tail hallucinations never
     pollute reported contact counts (while keeping continuous w_sb/w_hb/w_hp intact).
  2. Adds an 'Ordered only (w >= 0.25)' toggle (default ON) to the Contact Inspector
     Table and prints wHP inside Step 4 of index.html for 100% pocket-calculator parity.
  3. Re-runs Test A (fold_vip1_bzip61.zip) and Test B (--generate-demo).
"""

import ast
import re
import subprocess
import sys
from pathlib import Path

idx_path = Path("index.html")
cli_path = Path("unippi_cli.py")

# ==============================================================================
# PART 1: PATCH index.html
# ==============================================================================
html = idx_path.read_text(encoding="utf-8")

# 1. Gate displayed integer SB/HB/HP/pi/SS/clash counts by wij >= 0.25
old_js_cnt = "          nSb += effSb; nHb += effHb; nHp += effHp; nPi += pPi; nSs += pSs; nClash += effClash;"
new_js_cnt = "          if (wij >= 0.25) { nSb += effSb; nHb += effHb; nHp += effHp; nPi += pPi; nSs += pSs; nClash += effClash; }"
assert old_js_cnt in html, "Could not locate integer contact count line in index.html"
html = html.replace(old_js_cnt, new_js_cnt, 1)

# 2. Add wHP to Step 4 dG_phys derivation printout
old_step4_hp = "IC<sub>AA</sub>=${(ic.APOLAR_APOLAR||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å², wSB=${r.wSb||0}, wHB=${r.wHb||0}"
new_step4_hp = "IC<sub>AA</sub>=${(ic.APOLAR_APOLAR||0).toFixed(1)}, wHP=${(r.wHp||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å², wSB=${(r.wSb||0).toFixed(1)}, wHB=${(r.wHb||0).toFixed(1)}"
assert old_step4_hp in html, "Could not locate Step 4 dG_phys line in index.html"
html = html.replace(old_step4_hp, new_step4_hp, 1)

# 3. Add 'Ordered only (w >= 0.25)' toggle to Contact Inspector Table (checked by default)
old_filter_input = '<input type="text" id="contactFilterInput" placeholder="Filter (e.g. GLU, ARG, SB)..." style="padding:6px 10px; font-size:12px; border:1px solid var(--border-strong); border-radius:6px; font-family:var(--font-mono); width:185px; background:#fff;" oninput="renderContactInspectorTable(ALL_RESULTS[ACTIVE_INDEX])">'
new_filter_wrap = """<div style="display:flex; align-items:center; gap:8px; flex-wrap:wrap;">
              <label style="font-size:11.5px; font-weight:600; display:flex; align-items:center; gap:4px; cursor:pointer; color:var(--text-muted);">
                <input type="checkbox" id="hideLowConfContactsCheck" checked onchange="renderContactInspectorTable(ALL_RESULTS[ACTIVE_INDEX])">
                Ordered only (w ≥ 0.25)
              </label>
              """ + old_filter_input + "\n            </div>"
assert old_filter_input in html, "Could not locate contactFilterInput in index.html"
html = html.replace(old_filter_input, new_filter_wrap, 1)

old_filt_js = """      const q = (document.getElementById("contactFilterInput")?.value || "").trim().toUpperCase();
      const allC = (r.contacts || []).slice().sort((a,b) => a.minDist - b.minDist);
      const filt = allC.filter(c => !q || c.resA.toUpperCase().includes(q) || c.resB.toUpperCase().includes(q) || (q==="SB" && c.sb>0) || (q==="HB" && c.hb>0) || (q==="HP" && c.hp>0));"""
new_filt_js = """      const q = (document.getElementById("contactFilterInput")?.value || "").trim().toUpperCase();
      const hideLow = document.getElementById("hideLowConfContactsCheck")?.checked ?? true;
      const allC = (r.contacts || []).slice().sort((a,b) => a.minDist - b.minDist);
      const filt = allC.filter(c => (!hideLow || c.weight >= 0.25) && (!q || c.resA.toUpperCase().includes(q) || c.resB.toUpperCase().includes(q) || (q==="SB" && c.sb>0) || (q==="HB" && c.hb>0) || (q==="HP" && c.hp>0)));"""
assert old_filt_js in html, "Could not locate renderContactInspectorTable filter logic in index.html"
html = html.replace(old_filt_js, new_filt_js, 1)

idx_path.write_text(html, encoding="utf-8")
print("[PASS] Patched index.html cleanly.")

# ==============================================================================
# PART 2: DYNAMICALLY GATE DISPLAYED INTEGER COUNTS IN unippi_cli.py
# ==============================================================================
cli_lines = cli_path.read_text(encoding="utf-8").splitlines()
# Locate the continuous weighted sum line 'w_sb +=' and detect the weight variable name & integer counter lines directly above it
idx_wsb = next((i for i, l in enumerate(cli_lines) if "w_sb +=" in l), -1)
assert idx_wsb > 0, "Could not locate 'w_sb +=' in unippi_cli.py"

m_wvar = re.search(r'w_sb\s*\+=\s*([A-Za-z_][A-Za-z0-9_]*)\s*\*', cli_lines[idx_wsb])
assert m_wvar, f"Could not parse weight variable from: {cli_lines[idx_wsb]}"
w_var = m_wvar.group(1)

# Collect contiguous integer counter lines (n_sb, n_hb, n_hp, n_pi, n_ss, n_clash) immediately above idx_wsb
cnt_start = idx_wsb
while cnt_start > 0 and re.match(r'^\s*n_(?:sb|hb|hp|pi|ss|clash|salt|hbond|hydro)', cli_lines[cnt_start - 1]):
    cnt_start -= 1

assert cnt_start < idx_wsb, f"No integer counter lines found immediately above L{idx_wsb+1}: {cli_lines[idx_wsb-3:idx_wsb+2]}"

if f"if {w_var} >= 0.25:" not in cli_lines[cnt_start - 1]:
    ind = cli_lines[cnt_start][:len(cli_lines[cnt_start]) - len(cli_lines[cnt_start].lstrip())]
    gated_lines = [f"{ind}if {w_var} >= 0.25:"] + [f"    {l}" for l in cli_lines[cnt_start:idx_wsb]]
    cli_lines = cli_lines[:cnt_start] + gated_lines + cli_lines[idx_wsb:]
    cli_new = "\n".join(cli_lines) + "\n"
    ast.parse(cli_new)
    cli_path.write_text(cli_new, encoding="utf-8")
    print(f"[PASS] Gated {idx_wsb - cnt_start} integer contact counter line(s) by 'if {w_var} >= 0.25:' in unippi_cli.py (ast.parse OK).\n")
else:
    print("[INFO] Integer contact counters already gated in unippi_cli.py.\n")

# ==============================================================================
# PART 3: LIVE VERIFICATION (TEST A + TEST B)
# ==============================================================================
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
