#!/usr/bin/env python3
"""
Step 10:
  1. Gates resContrib in index.html by `if (wij >= 0.25)` (fixes Panel 3 A:LYS212 SB:0.4 leak).
  2. Formats UIS, LIS, and ipTM to 3 decimal places in #ensembleSummaryBanner (fixes +/- 0.00).
  3. Diagnoses the exact 0.04-0.11 kcal/mol CLI-vs-Browser delta on fold_vip1_bzip61_model_0
     by comparing unippi_cli.py's exact intermediates against the live browser output.
"""

import json
import re
from pathlib import Path

idx_path = Path("index.html")
cli_path = Path("unippi_cli.py")
html = idx_path.read_text(encoding="utf-8")

# 1. Fix Flaw #2: Gate resContrib by wij >= 0.25 so discarded tail contacts never color Panel 3
old_rc = """          for (const rk of [ra.key, rb.key]) {
            if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0 });
            const rc = resContrib.get(rk);
            rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*pPi;
          }"""
new_rc = """          if (wij >= 0.25) {
            for (const rk of [ra.key, rb.key]) {
              if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0 });
              const rc = resContrib.get(rk);
              rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*pPi;
            }
          }"""
assert old_rc in html, "Could not locate resContrib accumulation loop in index.html"
html = html.replace(old_rc, new_rc, 1)

# 2. Fix Flaw #3: 3-decimal Mean +/- SD for UIS, LIS & ipTM in #ensembleSummaryBanner
old_ens = """          const mSd = arr => {
            if (!arr.length) return "N/A";
            const m = arr.reduce((a,b)=>a+b,0)/arr.length;
            const sd = arr.length > 1 ? Math.sqrt(arr.reduce((a,b)=>a+Math.pow(b-m,2),0)/(arr.length-1)) : 0.0;
            return `${m.toFixed(2)} ± ${sd.toFixed(2)}`;
          };
          const fullPaeModels = ALL_RESULTS.filter(x => x.mode === "MODE_A_PRED_FULL");
          const lisSummary = fullPaeModels.length ? mSd(fullPaeModels.map(x=>x.LIS)) : "N/A (Coord Mode)";
          ensBanner.style.display = "block";
          ensBanner.innerHTML = `<strong>MULTI-MODEL ENSEMBLE SUMMARY (N = ${ALL_RESULTS.length} Models, Mean ± SD):</strong> UIS = <strong>${mSd(ALL_RESULTS.map(x=>x.UIS))}</strong> · ΔG<sub>uni</sub> = <strong>${mSd(ALL_RESULTS.map(x=>x.dG_uni_kcal))} kcal/mol</strong> · Local LIS = <strong>${lisSummary}</strong> · wBSA = <strong>${mSd(ALL_RESULTS.map(x=>x.wBSA_A2))} Å²</strong> · Confident Contacts = <strong>${mSd(ALL_RESULTS.map(x=>x.confContactCount ?? (x.contacts||[]).length))}</strong>`;"""
new_ens = """          const mSd = (arr, dec=2) => {
            if (!arr.length) return "N/A";
            const m = arr.reduce((a,b)=>a+b,0)/arr.length;
            const sd = arr.length > 1 ? Math.sqrt(arr.reduce((a,b)=>a+Math.pow(b-m,2),0)/(arr.length-1)) : 0.0;
            return `${m.toFixed(dec)} ± ${sd.toFixed(dec)}`;
          };
          const fullPaeModels = ALL_RESULTS.filter(x => x.mode === "MODE_A_PRED_FULL");
          const iptmSummary = fullPaeModels.length ? mSd(fullPaeModels.map(x=>x.ipTM), 3) : "N/A";
          const lisSummary = fullPaeModels.length ? mSd(fullPaeModels.map(x=>x.LIS), 3) : "N/A (Coord Mode)";
          ensBanner.style.display = "block";
          ensBanner.innerHTML = `<strong>MULTI-MODEL ENSEMBLE SUMMARY (N = ${ALL_RESULTS.length} Models, Mean ± SD):</strong> UIS = <strong>${mSd(ALL_RESULTS.map(x=>x.UIS), 3)}</strong> · ΔG<sub>uni</sub> = <strong>${mSd(ALL_RESULTS.map(x=>x.dG_uni_kcal), 2)} kcal/mol</strong> · ipTM = <strong>${iptmSummary}</strong> · Local LIS = <strong>${lisSummary}</strong> · wBSA = <strong>${mSd(ALL_RESULTS.map(x=>x.wBSA_A2), 1)} Å²</strong> · Confident Contacts = <strong>${mSd(ALL_RESULTS.map(x=>x.confContactCount ?? (x.contacts||[]).length), 1)}</strong>`;"""
assert old_ens in html, "Could not locate ensembleSummaryBanner block in index.html"
html = html.replace(old_ens, new_ens, 1)

idx_path.write_text(html, encoding="utf-8")
print("[PASS] Applied Flaw #2 (resContrib wij>=0.25 gate) & Flaw #3 (3-decimal UIS/LIS/ipTM ensemble summary) to index.html.\n")

# 3. Diagnose exact CLI vs Browser difference on fold_vip1_bzip61_model_0
rep_path = Path("./af3_vip1_bzip61_verified/unippi_full_report.json")
if rep_path.exists():
    rep = json.loads(rep_path.read_text(encoding="utf-8"))
    m0 = next((r for r in rep if r.get("complex_id", "").endswith("model_0")), None)
    m2 = next((r for r in rep if r.get("complex_id", "").endswith("model_2")), None)
    print("="*82)
    print("PARITY DIAGNOSTIC: CLI (unippi_full_report.json) vs LIVE BROWSER (model_0)")
    print("="*82)
    browser_m0 = {
        "UIS": 0.5653, "ipTM": 0.380, "LIS": 0.501, "pDockQ_PAE": 0.728,
        "dG_prodigy_kcal": -8.49, "dG_phys_kcal": -12.67, "dG_uni_kcal": -10.86,
        "Kd_nM": 10.92, "BSA_A2": 5171.5, "wBSA_A2": 4957.5,
        "mean_intf_pLDDT": 90.8, "mean_intf_PAE_A": 1.76, "mean_w_intf": 0.959,
        "salt_bridges": 11, "h_bonds": 12, "hydrophobic_contacts": 66, "intf_residues": 91
    }
    if m0:
        for k, bv in browser_m0.items():
            cv = m0.get(k, "N/A")
            diff = f"{abs(cv - bv):.4f}" if isinstance(cv, (int, float)) else "-"
            flag = " <-- DIFF" if isinstance(cv, (int, float)) and abs(cv - bv) > 1e-3 else ""
            print(f"{k:<22} | Browser: {str(bv):>10} | CLI: {str(cv):>10} | Delta: {diff:>8}{flag}")
    if m2:
        print("\n[model_2 check] CLI dG_PRODIGY =", m2.get("dG_prodigy_kcal"),
              "| dG_phys =", m2.get("dG_phys_kcal"), "| dG_uni =", m2.get("dG_uni_kcal"),
              "| ipTM =", m2.get("ipTM"), "| pDockQ_PAE =", m2.get("pDockQ_PAE"), "| UIS =", m2.get("UIS"))
        print("[model_2 check] Browser dG_uni = -10.94 | ipTM = 0.37 | pDockQ_PAE = 0.728 | UIS = 0.565")

# Also print PRODIGY_CLASS & NIS lines from unippi_cli.py so we see any mismatch immediately
cli_src = cli_path.read_text(encoding="utf-8")
print("\n" + "="*82)
print("CLI PRODIGY_CLASS / NIS / pair_cls SNIPPETS:")
print("="*82)
for i, line in enumerate(cli_src.splitlines(), 1):
    if any(tok in line for tok in ("PRODIGY_CLASS", "CHARGED_APOLAR", "POLAR_APOLAR", "nis_tot", "pct_nis", "dg_prodigy =", "p_pae", "pair_cls")):
        print(f"L{i:4d}: {line}")
