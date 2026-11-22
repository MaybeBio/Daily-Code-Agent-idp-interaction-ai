#!/usr/bin/env python3
"""
Step 9 (Final Phase 1 Edge-Case Seal before Browser .ZIP Drop-Test & Git Push):
  1. Hoists activeCutoff DOM lookup out of the inner residue-pair loop (faster 5-model .zip runs).
  2. Enables Panel 3 SASA burial rendering for NUCLEIC (DNA/RNA) chains (e.g. bZIP+DNA complexes).
  3. Uses an honest Coordinate-Mode pattern description when no PAE/ipTM JSON is provided (MODE_B/C).
  4. Preserves .pdb vs .cif extension inside downloadActiveRawAF3Pair().
  5. Audits on-disk index.html for live HTTPS CDN links & zero file:/// links + prints git status.
"""

import subprocess
from pathlib import Path

idx_path = Path("index.html")
html = idx_path.read_text(encoding="utf-8")

# 1. Hoist activeCutoff out of the inner contact loop
old_inner_cut = "          const activeCutoff = clamp(parseFloat(document.getElementById('cutoffInput')?.value || '5.5') || 5.5, 3.5, 8.0);\n"
assert old_inner_cut in html, "Could not locate inner activeCutoff line"
html = html.replace(old_inner_cut, "", 1)

anchor_top = "      const sha256_id = fnv1aHashHex([seqA, seqB].sort().join(\"::\"));"
assert anchor_top in html, "Could not locate sha256_id anchor in analyzePairJS()"
html = html.replace(
    anchor_top,
    anchor_top + "\n      const activeCutoff = clamp(parseFloat(document.getElementById('cutoffInput')?.value || '5.5') || 5.5, 3.5, 8.0);",
    1
)

# 2. Support NUCLEIC (DNA/RNA) chains in mutCatalog & Panel 3 labels
old_prot_only = '          if (r.entityType !== "PROTEIN") return;'
new_prot_nuc  = '          if (r.entityType !== "PROTEIN" && r.entityType !== "NUCLEIC") return;'
assert old_prot_only in html, "Could not locate PROTEIN-only mutCatalog guard"
html = html.replace(old_prot_only, new_prot_nuc, 1)

html = html.replace(
    '${THREE_TO_ONE[m.wtName]||"X"}',
    '${THREE_TO_ONE[m.wtName] || m.wtName.slice(-1) || "X"}'
)

# 3. Honest Coordinate-Mode pattern text (prevents calling X-ray PDBs 'IDR-Diluted low ipTM')
old_pat_end = """      } else if (wBsa >= 250.0 && dgUni > -5.0) {
        patId = "Pattern C"; patDesc = "Weak / Sub-Threshold Packing (Geometric burial without sufficient stabilizing contacts)";
      }"""
new_pat_end = old_pat_end + """
      if (mode !== "MODE_A_PRED_FULL" && (patId === "Pattern A" || patId === "Pattern B")) {
        patId = "Pattern A";
        patDesc = "Coordinate-Mode Supported Interactor (Interface packing pDockQ_coord and ΔG_uni support binding; no PAE JSON supplied)";
      }"""
assert old_pat_end in html, "Could not locate Pattern A-D closing branch"
html = html.replace(old_pat_end, new_pat_end, 1)

# 4. Preserve .pdb vs .cif extension in downloadActiveRawAF3Pair()
old_dl_ext = "      if (r.rawStructText) downloadTextFile(`${r.complex_id}.cif`, r.rawStructText);"
new_dl_ext = '      if (r.rawStructText) downloadTextFile(`${r.complex_id}.${r.rawStructFormat || "cif"}`, r.rawStructText);'
assert old_dl_ext in html, "Could not locate downloadActiveRawAF3Pair .cif line"
html = html.replace(old_dl_ext, new_dl_ext, 1)

idx_path.write_text(html, encoding="utf-8")

# 5. On-disk release hygiene audit
disk_h = idx_path.read_text(encoding="utf-8")
cdn_3dmol = [l.strip() for l in disk_h.splitlines() if "3Dmol-min.js" in l][0]
cdn_jszip = [l.strip() for l in disk_h.splitlines() if "jszip.min.js" in l][0]
assert "https://" in cdn_3dmol, f"3Dmol tag not HTTPS on disk: {cdn_3dmol}"
assert "https://" in cdn_jszip, f"JSZip tag not HTTPS on disk: {cdn_jszip}"
assert "file:///" not in disk_h, "Found stray file:/// URL on disk in index.html!"
assert "chuma" not in disk_h.lower(), "Found stray 'chuma' on disk in index.html!"

print("[PASS] Step 9 edge-case seal & on-disk release hygiene audit PASSED!")
print(f"  3Dmol CDN : {cdn_3dmol}")
print(f"  JSZip CDN : {cdn_jszip}")
print("\n[GIT STATUS]")
subprocess.run(["git", "status", "-s"], check=False)
