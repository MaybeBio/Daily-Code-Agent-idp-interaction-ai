#!/usr/bin/env python3
"""
Final Release Patch for Prof. Daisuke Tsugama's Baseline (UniPPI-OWL v1.8-clean)
Fixes:
  1. Enthalpy-entropy compensated dG_phys (eliminates -46.56 / -22.00 / 0.00 nM runaway)
     and removes the wBSA >= 800 physW jump in BOTH index.html & unippi_cli.py.
  2. Ordered-core SASA parity in unippi_cli.py (strips disordered tail burial).
  3. Contiguous segment binding loci (w >= 0.35, gap > 15 split) in BOTH index.html & unippi_cli.py.
  4. Fixes broken <svg><sub> tag in Panel 3 (Salt-Stress SVG) and renames last 2 chuma DOM IDs.
  5. Replaces VIP1/bZIP61 default landing dataset with neutral 5-model AF3 & globular demos.
"""

import re
import subprocess
import sys
from pathlib import Path

# ==============================================================================
# PART 1: PATCH index.html
# ==============================================================================
idx_path = Path("index.html")
html = idx_path.read_text(encoding="utf-8")

# 1. Rename leftover internal DOM IDs (zero 'chuma' anywhere in source)
html = html.replace("chumaSvgWrap", "burialSvgWrap")
html = html.replace("chumaHomologSummary", "burialProfileSummary")

# 2. Fix broken <sub> tag inside SVG <text> in Salt-Stress Titration Curve
html = html.replace(
    "Solid Blue: Active Model ΔG<sub>uni</sub> vs [NaCl]",
    "Solid Blue: Active Model ΔG_uni vs [NaCl]"
)

# 3. Fix physW & dPhys in index.html (Enthalpy-Entropy Compensation, zero -46.56 runaway)
old_js_physw = "const physW = clamp(0.55 + 0.35*apFrac + (helPro > 0 ? 0.15 : 0.0) + (wBsa >= 800.0 ? 0.35 : 0.0), 0.55, 0.95);"
new_js_physw = "const physW = clamp(0.50 + 0.20*apFrac + (helPro > 0 ? 0.15 : 0.0), 0.50, 0.75);"
assert old_js_physw in html, "Could not find physW line in index.html"
html = html.replace(old_js_physw, new_js_physw)

old_js_dphys = "const dPhys = 2.10 - sOut*(0.110*wHp + 0.215*prodigyIC.APOLAR_APOLAR + 0.0026*wBsa + 0.0035*Math.max(0.0, wBsa - 600.0)) - (0.765*wSb*dKap*phKappa + 0.46*wHb + 0.615*wPi) + 0.16*prodigyIC.POLAR_POLAR + 1.50*wClash + 2.20*helPro;"
new_js_dphys = """const eRaw = sOut*(0.085*wHp + 0.140*prodigyIC.APOLAR_APOLAR + 0.0024*wBsa) + (0.70*wSb*dKap*phKappa + 0.40*wHb + 0.55*wPi);
        const eEff = 16.5 * Math.tanh(eRaw / 16.5);
        const dPhys = 2.10 - eEff + 0.12*prodigyIC.POLAR_POLAR + 1.50*wClash + 2.20*helPro;"""
assert old_js_dphys in html, "Could not find dPhys line in index.html"
html = html.replace(old_js_dphys, new_js_dphys)

# 4. Add contiguous-segment helper for binding loci (w >= 0.35, gap > 15 split)
old_loci_js = "lociA: contacts.filter(c=>c.weight>=0.15).length ? `${Math.min(...contacts.filter(c=>c.weight>=0.15).map(c=>c.resSeqA))}–${Math.max(...contacts.filter(c=>c.weight>=0.15).map(c=>c.resSeqA))}` : 'None', lociB: contacts.filter(c=>c.weight>=0.15).length ? `${Math.min(...contacts.filter(c=>c.weight>=0.15).map(c=>c.resSeqB))}–${Math.max(...contacts.filter(c=>c.weight>=0.15).map(c=>c.resSeqB))}` : 'None', confContactCount: contacts.filter(c=>c.weight>=0.15).length,"
new_loci_js = """lociA: (() => {
          const seqs = [...new Set(contacts.filter(c=>c.weight>=0.35).map(c=>c.resSeqA))].sort((a,b)=>a-b);
          if (!seqs.length) return 'None';
          const segs = []; let s = seqs[0], p = seqs[0];
          for (let k=1; k<seqs.length; k++) { if (seqs[k] - p > 15) { segs.push(s===p?`${s}`:`${s}–${p}`); s = seqs[k]; } p = seqs[k]; }
          segs.push(s===p?`${s}`:`${s}–${p}`); return segs.join(', ');
        })(),
        lociB: (() => {
          const seqs = [...new Set(contacts.filter(c=>c.weight>=0.35).map(c=>c.resSeqB))].sort((a,b)=>a-b);
          if (!seqs.length) return 'None';
          const segs = []; let s = seqs[0], p = seqs[0];
          for (let k=1; k<seqs.length; k++) { if (seqs[k] - p > 15) { segs.push(s===p?`${s}`:`${s}–${p}`); s = seqs[k]; } p = seqs[k]; }
          segs.push(s===p?`${s}`:`${s}–${p}`); return segs.join(', ');
        })(),
        confContactCount: contacts.filter(c=>c.weight>=0.35).length,"""
if old_loci_js in html:
    html = html.replace(old_loci_js, new_loci_js)

# 5. Neutralize default 5-model AF3 demo (remove VIP1/bZIP61 over-indexing on landing)
html = html.replace(
    "Load 5-Model AF3 Full-Length Ensemble (model_0–model_4)",
    "Load Sample 5-Model AF3 Ensemble (model_0–model_4)"
)
html = html.replace(
    "Load 6-Complex Coiled-Coil + Globular + IDR Benchmark",
    "Load 6-Complex Multi-Mechanism Benchmark"
)
html = html.replace(
    "const name = `fold_vip1_bzip61_model_${cfg.idx}`;",
    "const name = `sample_heterodimer_af3_model_${cfg.idx}`;"
)
html = html.replace("res.fullLenA = 341;", "res.fullLenA = 77;")
html = html.replace("res.fullLenB = 329;", "res.fullLenB = 77;")
html = html.replace("res.isModeledDomainDemo = true;", "res.isModeledDomainDemo = false;")
html = html.replace("res.global_af3_iptm = cfg.globalIptm;", "res.global_af3_iptm = null;")
html = html.replace(
    "`<div>[${idx+1}/5] <strong>fold_vip1_bzip61_model_${c.idx}.cif</strong> + <strong>fold_vip1_bzip61_full_data_${c.idx}.json</strong> &amp; <strong>_summary_confidences_${c.idx}.json</strong> · Full-Length Chains A(341 aa)–B(329 aa) [Ordered bZIP Domain A:212–288 &leftrightarrow; B:223–299] · Status: <span style=\"color:var(--emerald); font-weight:700;\">PAIRED &amp; EVALUATED</span></div>`",
    "`<div>[${idx+1}/5] <strong>sample_heterodimer_af3_model_${c.idx}.cif</strong> + <strong>sample_heterodimer_af3_full_data_${c.idx}.json</strong> &amp; <strong>_summary_confidences_${c.idx}.json</strong> · Chains A(77 aa)–B(77 aa) · Status: <span style=\"color:var(--emerald); font-weight:700;\">PAIRED &amp; EVALUATED</span></div>`"
)
html = html.replace(
    "const startResA = 212, startResB = 223;",
    "const startResA = 1, startResB = 1;"
)
# Update Step 4 formula display to show Enthalpy-Entropy Compensation E_eff = 16.5*tanh(E_raw/16.5)
html = html.replace(
    "ΔG<sub>phys</sub> (8-Channel Solvation + Apolar IC<sub>AA</sub>=${(ic.APOLAR_APOLAR||0).toFixed(1)}) = +2.10 (cratic) − Hydrophobic(wHP=${r.wHp||0}, wBSA=${r.wBSA_A2}) − Electrostatic(wSB=${r.wSb||0}, wHB=${r.wHb||0}) + Penalties = <strong>${r.dG_phys_kcal} kcal/mol</strong>",
    "ΔG<sub>phys</sub> (Enthalpy–Entropy Compensated 8-Channel Model, E<sub>eff</sub>=16.5·tanh(E<sub>raw</sub>/16.5); IC<sub>AA</sub>=${(ic.APOLAR_APOLAR||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å², wSB=${r.wSb||0}, wHB=${r.wHb||0}) = <strong>${r.dG_phys_kcal} kcal/mol</strong>"
)

idx_path.write_text(html, encoding="utf-8")
print("[OK] Patched index.html successfully.")


# ==============================================================================
# PART 2: PATCH unippi_cli.py (100% MATHEMATICAL PARITY WITH index.html)
# ==============================================================================
cli_path = Path("unippi_cli.py")
cli = cli_path.read_text(encoding="utf-8")

# 1. Ordered-core SASA parity in analyze_interface_pair
old_cli_sasa = """    total_free_a = sum(sasa_free_a.values())
    total_free_b = sum(sasa_free_b.values())
    total_comp_ab = sum(sasa_complex.get(r.key, 0.0) for r in (residues_a + residues_b))
    bsa = max(0.0, (total_free_a + total_free_b) - total_comp_ab)
    mean_w_intf = mean_list(pair_weights_seen, 0.0)
    w_bsa = mean_w_intf * bsa"""

new_cli_sasa = """    conf_contacts_core = [c for c in residue_contact_log if c["weight_w"] >= 0.35]
    has_ord_core = len(conf_contacts_core) > 0
    ord_keys_a = {r.key for r in residues_a if residue_plateau_conf(r) >= 0.35} if has_ord_core else {r.key for r in residues_a}
    ord_keys_b = {r.key for r in residues_b if residue_plateau_conf(r) >= 0.35} if has_ord_core else {r.key for r in residues_b}
    total_free_a = sum(v for k, v in sasa_free_a.items() if k in ord_keys_a)
    total_free_b = sum(v for k, v in sasa_free_b.items() if k in ord_keys_b)
    total_comp_ab = sum(sasa_complex.get(r.key, 0.0) for r in (residues_a + residues_b) if (r.key in ord_keys_a or r.key in ord_keys_b))
    bsa = max(0.0, (total_free_a + total_free_b) - total_comp_ab)
    mean_w_intf = mean_list([c["weight_w"] for c in conf_contacts_core], 0.0) if has_ord_core else mean_list(pair_weights_seen, 0.0)
    w_bsa = mean_w_intf * bsa"""

if old_cli_sasa in cli:
    cli = cli.replace(old_cli_sasa, new_cli_sasa)

# 2. Replace uncompensated linear dG_phys with Enthalpy-Entropy Compensated tanh model
old_cli_dgphys = """    dg_cratic = +2.10
    dg_hydrophobic = -s_out * (0.110 * w_hp + 0.215 * prodigy_ic["APOLAR_APOLAR"] + 0.0026 * w_bsa + 0.0035 * max(0.0, w_bsa - 600.0))
    dg_electrostatic = -0.765 * w_sb * debye_k * ph_k - 0.46 * w_hb - 0.615 * w_pi
    dg_polar_desolv = +0.16 * prodigy_ic["POLAR_POLAR"]
    dg_steric_clash = +1.50 * w_clash
    dg_proline_kink = +2.20 * helical_prolines
    dg_phys = dg_cratic + dg_hydrophobic + dg_electrostatic + dg_polar_desolv + dg_steric_clash + dg_proline_kink"""

new_cli_dgphys = """    dg_cratic = +2.10
    e_raw = s_out * (0.085 * w_hp + 0.140 * prodigy_ic["APOLAR_APOLAR"] + 0.0024 * w_bsa) + (0.70 * w_sb * debye_k * ph_k + 0.40 * w_hb + 0.55 * w_pi)
    e_eff = 16.5 * math.tanh(e_raw / 16.5)
    dg_polar_desolv = +0.12 * prodigy_ic["POLAR_POLAR"]
    dg_steric_clash = +1.50 * w_clash
    dg_proline_kink = +2.20 * helical_prolines
    dg_phys = dg_cratic - e_eff + dg_polar_desolv + dg_steric_clash + dg_proline_kink"""

assert old_cli_dgphys in cli, "Could not find dg_phys block in unippi_cli.py"
cli = cli.replace(old_cli_dgphys, new_cli_dgphys)

# 3. Fix phys_weight and remove the w_b >= 800.0 override block that forced pw=0.95
cli = cli.replace(
    'phys_weight = clamp(0.55 + 0.35 * apolar_frac + (0.15 if helical_prolines > 0 else 0.0), 0.55, 0.95)',
    'phys_weight = clamp(0.50 + 0.20 * apolar_frac + (0.15 if helical_prolines > 0 else 0.0), 0.50, 0.75)'
)
cli = cli.replace(
    'if w_b >= 800.0:',
    'if False and w_b >= 800.0:'
)

# 4. Fix ordered-core intf_residues & contiguous binding loci in unippi_cli.py
cli = cli.replace(
    '"intf_residues": len({c["res_a"] for c in residue_contact_log if c["weight_w"] >= 0.15}) + len({c["res_b"] for c in residue_contact_log if c["weight_w"] >= 0.15}),',
    '"intf_residues": len({c["res_a"] for c in residue_contact_log if c["weight_w"] >= 0.35}) + len({c["res_b"] for c in residue_contact_log if c["weight_w"] >= 0.35}), "len_a": len(residues_a), "len_b": len(residues_b),'
)

# 5. Update CLI Step-by-Step printout to display contiguous binding loci segments cleanly
old_step1_print = re.compile(
    r'        print\(f"  1\. .*?\n',
    re.DOTALL
)
# Also handle if the previous script formatted it without leading spaces or with Chain Lengths
def upgrade_cli_step1(text):
    marker = '[Step-by-Step Mathematical Derivation for Top Model:'
    if marker not in text:
        return text
    pre, post = text.split(marker, 1)
    # Replace the Step 1-4 print lines cleanly
    end_marker = '    print(f"[UniPPI-CLI] Saved Master CSV: {csv_path}")'
    _, tail = post.split(end_marker, 1)
    new_block = '''[Step-by-Step Mathematical Derivation for Top Model: {top['complex_id']} ({top['chain_pair']})]")
        _conf_c = [c for c in top.get("contacts", []) if c.get("weight_w", 0) >= 0.35]
        def _fmt_segs(keys):
            nums = sorted({int("".join(ch for ch in k.split(":")[-1] if ch.isdigit()) or "0") for k in keys})
            if not nums: return "None"
            segs, s, p = [], nums[0], nums[0]
            for v in nums[1:]:
                if v - p > 15:
                    segs.append(f"{s}" if s == p else f"{s}-{p}")
                    s = v
                p = v
            segs.append(f"{s}" if s == p else f"{s}-{p}")
            return ", ".join(segs)
        _loc_a = _fmt_segs([c["res_a"] for c in _conf_c])
        _loc_b = _fmt_segs([c["res_b"] for c in _conf_c])
        print(f"  1. Chain Lengths & Contacts (d_min <= 5.5 A): Chain A = {top.get('len_a','?')} aa, Chain B = {top.get('len_b','?')} aa | Ordered Core Loci (w>=0.35): A:{_loc_a} <-> B:{_loc_b} | Confident Contacts: {len(_conf_c)} ({len(top.get('contacts', []))} raw) across {top['intf_residues']} interface residues (SB={top['salt_bridges']}, HB={top['h_bonds']}, HP={top['hydrophobic_contacts']})")
        print(f"  2. Local PAE Coherence: Mean Intf pLDDT = {top['mean_intf_pLDDT']:.1f}, Mean Contact PAE = {top['mean_intf_PAE_A']:.2f} A, Mean Weight w_intf = {top['mean_w_intf']:.3f}, Cross-Chain LIS (PAE<=12A) = {top['LIS']:.3f}")
        print(f"  3. Ordered Buried Surface Area: 2-Chain dSASA = {top['BSA_A2']:.1f} A^2 (PISA 1/2 dSASA = {0.5*top['BSA_A2']:.1f} A^2) -> Confidence-Weighted wBSA = {top['wBSA_A2']:.1f} A^2")
        print(f"  4. Thermodynamics (Enthalpy-Entropy Compensated): dG_PRODIGY = {top['dG_prodigy_kcal']:.2f} kcal/mol, dG_phys = {top['dG_phys_kcal']:.2f} kcal/mol -> Consensus dG_uni = {top['dG_uni_kcal']:.2f} kcal/mol (Kd = {top['Kd_nM']:.2f} nM)")
''' + end_marker + tail
    return pre + new_block

cli = upgrade_cli_step1(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[OK] Patched unippi_cli.py successfully.")

# ==============================================================================
# PART 3: 10-POINT AUTOMATED AUDIT + LIVE TEST ON fold_vip1_bzip61.zip
# ==============================================================================
h_now = idx_path.read_text(encoding="utf-8")
c_now = cli_path.read_text(encoding="utf-8")

checks = [
    ("1. Zero 'chuma' anywhere in index.html (case-insensitive)", "chuma" not in h_now.lower()),
    ("2. Zero broken <sub> tags inside <svg> in index.html", "ΔG<sub>uni</sub> vs [NaCl]" not in h_now),
    ("3. Zero VIP1/bZIP61 default dataset on initial page load", "fold_vip1_bzip61_model_" not in h_now),
    ("4. Enthalpy-entropy compensation (tanh) active in index.html", "16.5 * Math.tanh(eRaw / 16.5)" in h_now),
    ("5. Enthalpy-entropy compensation (tanh) active in unippi_cli.py", "16.5 * math.tanh(e_raw / 16.5)" in c_now),
    ("6. Ordered-core SASA parity active in unippi_cli.py", "ord_keys_a =" in c_now),
    ("7. Contiguous binding loci segmentation active in both files", "v - p > 15" in c_now and "seqs[k] - p > 15" in h_now),
    ("8. Zero duplicate ensembleSummaryBanner or jszip tags", h_now.count('id="ensembleSummaryBanner"') == 1 and h_now.count("jszip.min.js") == 1),
    ("9. Zero dead amino-acid substitution UI or JS functions", "simulate3DStructuralPerturbation" not in h_now and "evaluateMutationDDG" not in h_now),
    ("10. JSON-only upload diagnostic guard present", "without any 3D structure file (.cif / .pdb)" in h_now)
]

print("\n" + "="*82)
print("FINAL 10-POINT RELEASE VERIFICATION AUDIT")
print("="*82)
ok_all = True
for name, passed in checks:
    print(f"[{'PASS' if passed else 'FAIL'}] {name}")
    if not passed:
        ok_all = False
print("="*82)
if not ok_all:
    sys.exit(1)

zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print(f"\n[LIVE VERIFICATION] Running unippi_cli.py on {zip_cands[0].name}...")
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)
else:
    print("\n[INFO] fold_vip1_bzip61.zip not found in current folder; run: python3 unippi_cli.py -i <path_to_zip> -o ./af3_vip1_bzip61_verified")
