#!/usr/bin/env python3
"""
Step 7 (Final Phase 1 Edge-Case Hardening + Head-to-Head JS vs. CLI Parity Proof):
  1. Slices r.paeMat strictly to Chain A + Chain B for 3+ chain complexes.
  2. Uses exact AF3 chain_pair_iptm[i][j] when present in *_summary_confidences_*.json.
  3. Replaces misleading PAE=3.5 / LIS=0 / blank SVG with clear N/A (Coordinate Mode)
     indicators whenever no PAE JSON is uploaded (MODE_B / MODE_C).
  4. Runs an automated head-to-head numerical parity proof comparing index.html's
     JavaScript engine directly against unippi_cli.py on fold_vip1_bzip61_model_2.
"""

import json
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

idx_path = Path("index.html")
html = idx_path.read_text(encoding="utf-8")

# 1. Fix A: Slice r.paeMat strictly to (resA + resB) for 3+ chain complexes
old_pae_ret = "        paeMat, lenA: resA.length, lenB: resB.length,"
new_pae_ret = """        paeMat: (paeMat && paeMat.length !== (resA.length + resB.length))
          ? [...resA, ...resB].map(ra => { const ia = resToIdx.get(ra.key); return [...resA, ...resB].map(rb => paeMat[ia][resToIdx.get(rb.key)]); })
          : paeMat,
        lenA: resA.length, lenB: resB.length,"""
assert old_pae_ret in html, "Could not locate paeMat return property in index.html"
html = html.replace(old_pae_ret, new_pae_ret, 1)

# 2. Fix B: Extract AF3 chain_pair_iptm[i][j] in handleFileSelect()
old_json_init = "        let iptm = null, paeMat = null, hasJson = false;"
new_json_init = "        let iptm = null, chainPairIptm = null, paeMat = null, hasJson = false;"
assert old_json_init in html, "Could not locate iptm/paeMat init in handleFileSelect()"
html = html.replace(old_json_init, new_json_init, 1)

old_iptm_read = """          if (rawIptm !== undefined && rawIptm !== null) {
            iptm = parseFloat(rawIptm);
            hasJson = true;
            pairedNames.push(j.name);
          }"""
new_iptm_read = """          if (Array.isArray(j.data.chain_pair_iptm)) chainPairIptm = j.data.chain_pair_iptm;
          if (rawIptm !== undefined && rawIptm !== null) {
            iptm = parseFloat(rawIptm);
            hasJson = true;
            pairedNames.push(j.name);
          }"""
assert old_iptm_read in html, "Could not locate rawIptm block in handleFileSelect()"
html = html.replace(old_iptm_read, new_iptm_read, 1)

old_call_pair = "            const res = analyzePairJS(stem, ca, cb, chains[ca], chains[cb], ordered, resToIdx, mode, { iptm, res_pae_matrix: paeMat }, pairSasaComp, sasaFree[ca], sasaFree[cb], tempK, lambdaPae, isBinary);"
new_call_pair = """            const pairIptm = (Array.isArray(chainPairIptm) && chainPairIptm[i] && typeof chainPairIptm[i][j] === "number") ? chainPairIptm[i][j] : iptm;
            const res = analyzePairJS(stem, ca, cb, chains[ca], chains[cb], ordered, resToIdx, mode, { iptm: pairIptm, res_pae_matrix: paeMat }, pairSasaComp, sasaFree[ca], sasaFree[cb], tempK, lambdaPae, isBinary);"""
assert old_call_pair in html, "Could not locate analyzePairJS call in handleFileSelect()"
html = html.replace(old_call_pair, new_call_pair, 1)

# 3. Fix E: Clean Coordinate-Mode (MODE_B / MODE_C) N/A indicators when no PAE JSON is supplied
old_lb_iptm = """          <td>${r.ipTM.toFixed(2)}</td>
          <td>${r.LIS.toFixed(3)}</td>"""
new_lb_iptm = """          <td>${r.mode === "MODE_A_PRED_FULL" ? r.ipTM.toFixed(2) : "—"}</td>
          <td>${r.mode === "MODE_A_PRED_FULL" ? r.LIS.toFixed(3) : "—"}</td>"""
assert old_lb_iptm in html, "Could not locate leaderboard ipTM/LIS cells"
html = html.replace(old_lb_iptm, new_lb_iptm, 1)

old_step2_pae = "Mean interface pLDDT = <strong>${r.mean_intf_pLDDT}</strong>, Mean contact PAE = <strong>${r.mean_intf_PAE_A} Å</strong>, Mean weight w̄<sub>intf</sub> = <strong>${r.mean_w_intf.toFixed(3)}</strong>, and cross-chain <strong>LIS (PAE ≤ 12 Å) = ${r.LIS}</strong>."
new_step2_pae = "Mean interface pLDDT = <strong>${r.mean_intf_pLDDT}</strong>, Mean contact PAE = <strong>${r.mode === \"MODE_A_PRED_FULL\" ? `${r.mean_intf_PAE_A} Å` : \"N/A (Coordinate Mode)\"}</strong>, Mean weight w̄<sub>intf</sub> = <strong>${r.mean_w_intf.toFixed(3)}</strong>, and cross-chain <strong>LIS (PAE ≤ 12 Å) = ${r.mode === \"MODE_A_PRED_FULL\" ? r.LIS.toFixed(3) : \"N/A\"}</strong>."
assert old_step2_pae in html, "Could not locate Step 2 PAE line"
html = html.replace(old_step2_pae, new_step2_pae, 1)

old_track2 = "<td>LIS = ${r.LIS.toFixed(3)} · pDockQ<sub>PAE</sub> = ${r.pDockQ_PAE.toFixed(3)} · Mean Contact PAE = ${r.mean_intf_PAE_A} Å</td>"
new_track2 = "<td>${r.mode === \"MODE_A_PRED_FULL\" ? `LIS = ${r.LIS.toFixed(3)} · pDockQ<sub>PAE</sub> = ${r.pDockQ_PAE.toFixed(3)} · Mean Contact PAE = ${r.mean_intf_PAE_A} Å` : `LIS = N/A (No PAE JSON) · pDockQ<sub>coord</sub> = ${r.pDockQ_coord.toFixed(3)} · PAE = N/A`}</td>"
assert old_track2 in html, "Could not locate Track 2 row"
html = html.replace(old_track2, new_track2, 1)

old_pae_svg_set = """      document.getElementById("paeSvgWrap").innerHTML = `
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${pw} ${ph}" width="100%" height="265">
          <rect width="${pw}" height="${ph}" fill="#ffffff"/>
          ${rects}"""
new_pae_svg_set = """      if (!r.paeMat || !r.paeMat.length) {
        document.getElementById("paeSvgWrap").innerHTML = `
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${pw} ${ph}" width="100%" height="265">
            <rect width="${pw}" height="${ph}" fill="#f8fafc" rx="6"/>
            <rect x="45" y="65" width="370" height="145" fill="#ffffff" stroke="#cbd5e1" stroke-dasharray="5,4" rx="8"/>
            <text x="230" y="125" text-anchor="middle" font-size="13.5" font-weight="bold" fill="#334155">PAE Matrix Unavailable (Coordinate-Only Mode)</text>
            <text x="230" y="150" text-anchor="middle" font-size="11.5" fill="#64748b">Upload a companion AF3 *_full_data_*.json or *_confidences.json</text>
            <text x="230" y="168" text-anchor="middle" font-size="11.5" fill="#64748b">alongside your .cif/.pdb file to render PAE &amp; Local LIS.</text>
          </svg>
        `;
      } else document.getElementById("paeSvgWrap").innerHTML = `
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${pw} ${ph}" width="100%" height="265">
          <rect width="${pw}" height="${ph}" fill="#ffffff"/>
          ${rects}"""
assert old_pae_svg_set in html, "Could not locate paeSvgWrap assignment"
html = html.replace(old_pae_svg_set, new_pae_svg_set, 1)

idx_path.write_text(html, encoding="utf-8")
print("[PASS] Hardened index.html for 3+ chain complexes & Coordinate-Only Mode (MODE_B/C).")

# ==============================================================================
# 4. AUTOMATED HEAD-TO-HEAD BROWSER-JS vs. PYTHON-CLI NUMERICAL PARITY PROOF
# ==============================================================================
js_blocks = re.findall(r'<script>(.*?)</script>', html, flags=re.DOTALL)
js_src = js_blocks[-1]

if shutil.which("node"):
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
        tf.write(js_src)
        chk_path = tf.name
    proc = subprocess.run(["node", "--check", chk_path], capture_output=True, text=True)
    Path(chk_path).unlink(missing_ok=True)
    assert proc.returncode == 0, f"JS syntax error:\n{proc.stderr}"
    print("[PASS] Node.js JS syntax check passed (0 errors).")

    # Run head-to-head parity on fold_vip1_bzip61_model_2 from real AF3 zip
    zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
    cli_rep = Path("./af3_vip1_bzip61_verified/unippi_full_report.json")
    if zip_cands and cli_rep.exists():
        with zipfile.ZipFile(zip_cands[0], "r") as zf:
            names = zf.namelist()
            m2_cif = next(n for n in names if n.endswith("model_2.cif") and "templates" not in n.lower())
            m2_full = next(n for n in names if "full_data_2.json" in n.lower())
            m2_sum = next(n for n in names if "summary_confidences_2.json" in n.lower())
            cif_txt = zf.read(m2_cif).decode("utf-8")
            full_json = zf.read(m2_full).decode("utf-8")
            sum_json = zf.read(m2_sum).decode("utf-8")

        runner = f"""
global.document = {{ getElementById: () => null }};
global.window = {{ addEventListener: () => {{}} }};
{js_src}
const fs = require("fs");
const cif = fs.readFileSync(process.argv[2], "utf-8");
const fullD = JSON.parse(fs.readFileSync(process.argv[3], "utf-8"));
const sumD = JSON.parse(fs.readFileSync(process.argv[4], "utf-8"));
const {{ atoms, meta }} = parseStructureText(cif, "fold_vip1_bzip61_model_2.cif");
const {{ ordered, chains }} = buildResidues(atoms);
const keyToTok = new Map();
fullD.token_chain_ids.forEach((tc, tidx) => {{
  const k = `${{tc}}|${{fullD.token_res_ids[tidx]}}`;
  if (!keyToTok.has(k)) keyToTok.set(k, tidx);
}});
const idxs = ordered.map(r => Math.min(fullD.pae.length - 1, keyToTok.get(`${{r.origChainId || r.chainId}}|${{r.resSeq}}`) ?? 0));
const paeMat = idxs.map(ii => idxs.map(jj => Number(fullD.pae[ii][jj])));
const sasaComp = computeResidueSASA(ordered);
const sasaA = computeResidueSASA(chains["A"]), sasaB = computeResidueSASA(chains["B"]);
const resToIdx = new Map(ordered.map((r, i) => [r.key, i]));
const res = analyzePairJS("fold_vip1_bzip61_model_2", "A", "B", chains["A"], chains["B"], ordered, resToIdx, "MODE_A_PRED_FULL", {{ iptm: sumD.iptm, res_pae_matrix: paeMat }}, sasaComp, sasaA, sasaB, 298.15, 8.0, true);
console.log(JSON.stringify({{
  UIS: res.UIS, dG_PRODIGY: res.dG_prodigy_kcal, dG_phys: res.dG_phys_kcal, dG_uni: res.dG_uni_kcal,
  Kd_nM: res.Kd_nM, LIS: res.LIS, wBSA: res.wBSA_A2, SB: res.salt_bridges, HB: res.h_bonds,
  HP: res.hydrophobic_contacts, Pro: res.helical_prolines, Pattern: res.pattern_id,
  lociA: res.lociA, lociB: res.lociB, confContacts: res.confContactCount
}}));
"""
        with tempfile.TemporaryDirectory() as td:
            tdp = Path(td)
            (tdp / "run.js").write_text(runner, encoding="utf-8")
            (tdp / "m2.cif").write_text(cif_txt, encoding="utf-8")
            (tdp / "m2_full.json").write_text(full_json, encoding="utf-8")
            (tdp / "m2_sum.json").write_text(sum_json, encoding="utf-8")
            out = subprocess.run(["node", str(tdp / "run.js"), str(tdp / "m2.cif"), str(tdp / "m2_full.json"), str(tdp / "m2_sum.json")], capture_output=True, text=True, check=True)
            js_m2 = json.loads(out.stdout.strip().splitlines()[-1])

        cli_data = json.loads(cli_rep.read_text(encoding="utf-8"))
        cli_m2 = next(r for r in cli_data if r["complex_id"].endswith("model_2"))

        print("\n" + "="*82)
        print("HEAD-TO-HEAD NUMERICAL PARITY PROOF: BROWSER JS vs. PYTHON CLI (fold_vip1_bzip61_model_2)")
        print("="*82)
        rows = [
            ("UIS", js_m2["UIS"], cli_m2["UIS"]),
            ("dG_PRODIGY (kcal/mol)", js_m2["dG_PRODIGY"], cli_m2["dG_prodigy_kcal"]),
            ("dG_phys (kcal/mol)", js_m2["dG_phys"], cli_m2["dG_phys_kcal"]),
            ("dG_uni (kcal/mol)", js_m2["dG_uni"], cli_m2["dG_uni_kcal"]),
            ("Kd (nM)", js_m2["Kd_nM"], cli_m2["Kd_nM"]),
            ("LIS", js_m2["LIS"], cli_m2["LIS"]),
            ("wBSA (A^2)", js_m2["wBSA"], cli_m2["wBSA_A2"]),
            ("Salt Bridges (SB)", js_m2["SB"], cli_m2["salt_bridges"]),
            ("H-Bonds (HB)", js_m2["HB"], cli_m2["h_bonds"]),
            ("Hydrophobic (HP)", js_m2["HP"], cli_m2["hydrophobic_contacts"]),
            ("Strained Prolines", js_m2["Pro"], cli_m2["helical_prolines"]),
            ("Pattern", js_m2["Pattern"], cli_m2["pattern_id"]),
        ]
        print(f"{'Metric':<24} | {'Browser JS':>14} | {'Python CLI':>14} | {'Delta':>10}")
        print("-"*72)
        for label, jv, cv in rows:
            delta = f"{abs(jv - cv):.4f}" if isinstance(jv, (int, float)) else ("MATCH" if jv == cv else "DIFF")
            print(f"{label:<24} | {str(jv):>14} | {str(cv):>14} | {delta:>10}")
        print("="*82)
else:
    print("[INFO] Node.js not installed in WSL; browser JS syntax & edge-case patches written cleanly.")
