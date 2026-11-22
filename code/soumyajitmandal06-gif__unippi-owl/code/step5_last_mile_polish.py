#!/usr/bin/env python3
"""
Step 5 (Last-Mile Polish before Phase 1 Git Push & Prof. Daisuke Email):
  1. Fixes the < 50 aa PRODIGY key-sorting bug in index.html (unconditional remap).
  2. Adds live per-model progress banner + event-loop yield in handleFileSelect()
     so large 5-model AF3 .zip uploads never freeze the browser UI.
  3. Removes unsupported 'unippi_full_report.json' dropzone copy, auto-hides
     Download Ligand SDF when no ligand exists, and removes dead .lab-header CSS.
  4. Adds exact PRODIGY/physical blend weights to unippi_cli.py Step 4 printout.
  5. Verifies live HTTPS CDN links for 3Dmol & JSZip on disk (ready for GitHub Pages).
"""

import ast
import re
import subprocess
import sys
from pathlib import Path

idx_path = Path("index.html")
cli_path = Path("unippi_cli.py")
html = idx_path.read_text(encoding="utf-8")
cli = cli_path.read_text(encoding="utf-8")

# 1. Fix the < 50 aa PRODIGY key-sorting bug in index.html
old_prodigy_remap = """          if (Math.min(resA.length, resB.length) >= 50) {
            if (pairCls === "APOLAR_CHARGED") pairCls = "CHARGED_APOLAR";
            else if (pairCls === "APOLAR_POLAR") pairCls = "POLAR_APOLAR";
          }"""
new_prodigy_remap = """          if (pairCls === "APOLAR_CHARGED") pairCls = "CHARGED_APOLAR";
          else if (pairCls === "APOLAR_POLAR") pairCls = "POLAR_APOLAR";"""
assert old_prodigy_remap in html, "Could not locate <50 aa PRODIGY remap bug in index.html"
html = html.replace(old_prodigy_remap, new_prodigy_remap, 1)

# 2. Add live per-model progress indicator + frame yield in handleFileSelect()
old_loop_hdr = "      for (const sf of structEntries) {"
new_loop_hdr = """      const statusBan = document.getElementById("statusBanner");
      const statusTxt = document.getElementById("statusText");
      if (statusBan) statusBan.style.display = "flex";
      for (let sfIdx = 0; sfIdx < structEntries.length; sfIdx++) {
        const sf = structEntries[sfIdx];
        if (statusTxt) statusTxt.textContent = `Evaluating model ${sfIdx + 1}/${structEntries.length}: ${sf.name}...`;
        await new Promise(resolve => setTimeout(resolve, 15));"""
assert old_loop_hdr in html, "Could not locate structEntries loop in handleFileSelect()"
html = html.replace(old_loop_hdr, new_loop_hdr, 1)

# 3. Clean dropzone copy, context-aware Ligand SDF button, and dead .lab-header CSS
html = html.replace(" (or precomputed <code>unippi_full_report.json</code>)", "")
html = html.replace(
    '<button class="btn btn-secondary btn-sm" onclick="downloadLigandSDF()">Download Ligand SDF</button>',
    '<button id="btnDownloadLigandSDF" class="btn btn-secondary btn-sm" style="display:none;" onclick="downloadLigandSDF()">Download Ligand SDF</button>'
)
html = html.replace(
    "      renderContactInspectorTable(r);",
    """      const ligBtn = document.getElementById("btnDownloadLigandSDF");
      if (ligBtn) ligBtn.style.display = (r.sdfText && r.sdfText.length > 0) ? "inline-flex" : "none";
      renderContactInspectorTable(r);"""
)
html = re.sub(r'\n    \.lab-header \{.*?\n    \.lab-nav a:hover \{[^\n]+\}\n', '\n', html, flags=re.DOTALL, count=1)

# Ensure HTTPS CDN tags are intact on disk for GitHub Pages deployment
if "3Dmol-min.js" in html and "https://" not in html.split("3Dmol-min.js")[0].splitlines()[-1]:
    html = re.sub(r'<script src="[^"]*3Dmol-min\.js[^"]*"></script>', '<script src="https://3Dmol.org/build/3Dmol-min.js"></script>', html)
if "jszip.min.js" in html and "https://" not in html.split("jszip.min.js")[0].splitlines()[-1]:
    html = re.sub(r'<script src="[^"]*jszip\.min\.js[^"]*"></script>', '<script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>', html)

idx_path.write_text(html, encoding="utf-8")

# 4. Add exact blend weights into unippi_cli.py Step 4 printout (if present)
old_cli_step4 = "dG_PRODIGY = {top['dG_prodigy_kcal']:.2f} kcal/mol, dG_phys = {top['dG_phys_kcal']:.2f} kcal/mol -> Consensus dG_uni = {top['dG_uni_kcal']:.2f} kcal/mol"
new_cli_step4 = "dG_PRODIGY = {top['dG_prodigy_kcal']:.2f} kcal/mol, dG_phys = {top['dG_phys_kcal']:.2f} kcal/mol (blend: {round(100*(top['dG_uni_kcal']-top['dG_phys_kcal'])/(top['dG_prodigy_kcal']-top['dG_phys_kcal'])) if abs(top['dG_prodigy_kcal']-top['dG_phys_kcal'])>1e-3 else 0}% PRODIGY / {100 - (round(100*(top['dG_uni_kcal']-top['dG_phys_kcal'])/(top['dG_prodigy_kcal']-top['dG_phys_kcal'])) if abs(top['dG_prodigy_kcal']-top['dG_phys_kcal'])>1e-3 else 0)}% phys) -> Consensus dG_uni = {top['dG_uni_kcal']:.2f} kcal/mol"
if old_cli_step4 in cli:
    cli = cli.replace(old_cli_step4, new_cli_step4, 1)
    ast.parse(cli)
    cli_path.write_text(cli, encoding="utf-8")

# 5. Final Verification Assertions
h = idx_path.read_text(encoding="utf-8")
assert "Math.min(resA.length, resB.length) >= 50" not in h, "FAIL: <50 aa bug still present"
assert "unippi_full_report.json" not in h, "FAIL: unsupported dropzone copy still present"
assert ".lab-header {" not in h, "FAIL: dead .lab-header CSS still present"
assert 'id="btnDownloadLigandSDF"' in h, "FAIL: context-aware ligand button missing"
assert "https://" in [l for l in h.splitlines() if "3Dmol-min.js" in l][0], "FAIL: 3Dmol CDN not HTTPS"
assert "https://" in [l for l in h.splitlines() if "jszip.min.js" in l][0], "FAIL: JSZip CDN not HTTPS"
assert "Evaluating model ${sfIdx + 1}/${structEntries.length}" in h, "FAIL: async progress feedback missing"

print("[PASS] All 7 last-mile polish assertions PASSED!")

zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print(f"\n[FINAL CLI VERIFICATION] unippi_cli.py -i {zip_cands[0].name}")
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)
