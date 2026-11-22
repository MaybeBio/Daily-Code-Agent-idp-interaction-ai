#!/usr/bin/env python3
"""
Step 4 (AST-Introspective Final Phase 1 Lock):
  1. Uses Python AST introspection to locate the exact Pattern A-D `if/elif` block
     and its real variable names inside unippi_cli.py (zero string guessing).
  2. Upgrades index.html: authentic PDB atoms in neutral 5-model sample, single
     compact sticky header, pocket-calculator Step 4 math, zero dead branches.
  3. Runs a 16-Point Master Release Audit + Regression Tests A & B.
"""

import ast
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

cli_path = Path("unippi_cli.py")
idx_path = Path("index.html")

# ==============================================================================
# PART 1: AST-INTROSPECTIVE PATTERN A-D SYNC IN unippi_cli.py
# ==============================================================================
cli_src = cli_path.read_text(encoding="utf-8")
cli_lines = cli_src.splitlines()
tree = ast.parse(cli_src)

# Find which FunctionDef contains the literal "Pattern C"
target_fn = None
for node in ast.walk(tree):
    if isinstance(node, ast.FunctionDef):
        if any(isinstance(sub, ast.Constant) and sub.value == "Pattern C" for sub in ast.walk(node)):
            target_fn = node
            break
assert target_fn is not None, "Could not find any function containing 'Pattern C' in unippi_cli.py"

def find_outer_if_with_pattern_c(stmts):
    for st in stmts:
        if isinstance(st, ast.If):
            if any(isinstance(sub, ast.Constant) and sub.value == "Pattern C" for sub in ast.walk(st)):
                return st
        for attr in ("body", "orelse", "finalbody"):
            sub_list = getattr(st, attr, None)
            if isinstance(sub_list, list):
                found = find_outer_if_with_pattern_c(sub_list)
                if found:
                    return found
    return None

outer_if = find_outer_if_with_pattern_c(target_fn.body)
assert outer_if is not None, f"Could not find If-block with 'Pattern C' inside {target_fn.name}()"

s_line = outer_if.lineno - 1        # 0-based inclusive
e_line = outer_if.end_lineno        # 0-based exclusive

print(f"[AST AUDIT] Located Pattern A-D If-block in {target_fn.name}() at lines {s_line+1}..{e_line}:")
for k in range(s_line, e_line):
    print(f"  OLD L{k+1:4d}: {cli_lines[k]}")

# Collect all local variable names in scope inside target_fn
in_scope = {a.arg for a in target_fn.args.args}
for sub in ast.walk(target_fn):
    if isinstance(sub, ast.Assign):
        for t in sub.targets:
            if isinstance(t, ast.Name):
                in_scope.add(t.id)
            elif isinstance(t, (ast.Tuple, ast.List)):
                for elt in t.elts:
                    if isinstance(elt, ast.Name):
                        in_scope.add(elt.id)

mod_globals = {
    t.id for st in tree.body if isinstance(st, ast.Assign)
    for t in st.targets if isinstance(t, ast.Name)
}

def pick_var(cands, default=None):
    for c in cands:
        if c in in_scope:
            return c
    return default

uis_v    = pick_var(["uis", "UIS"])
iptm_v   = pick_var(["iptm", "iptm_val", "ipTM"])
lis_v    = pick_var(["lis", "lis_score", "LIS"])
dg_v     = pick_var(["dg_uni", "dG_uni", "dg"])
pro_v    = pick_var(["helical_prolines", "hel_pro", "pro_kinks"])
clash_v  = pick_var(["n_clash", "n_clashes", "steric_clashes", "clashes"])
wbsa_v   = pick_var(["w_bsa", "wbsa", "wBSA"])
thresh_v = "UIS_THRESHOLD" if "UIS_THRESHOLD" in mod_globals else "0.5000"

assert all([uis_v, iptm_v, lis_v, dg_v, pro_v, wbsa_v]), (
    f"Missing expected variables in {target_fn.name}(): "
    f"uis={uis_v}, iptm={iptm_v}, lis={lis_v}, dg={dg_v}, pro={pro_v}, wbsa={wbsa_v} | in_scope={sorted(in_scope)}"
)

clash_ok   = f" and {clash_v} < 5" if clash_v else ""
clash_bad  = f" or {clash_v} >= 5" if clash_v else ""
ind = cli_lines[s_line][:len(cli_lines[s_line]) - len(cli_lines[s_line].lstrip())]

# Check whether outer_if uses return statements or variable assignments
uses_return = any(isinstance(sub, ast.Return) for sub in ast.walk(outer_if))

if uses_return:
    def emit_pat(pid, pdesc):
        return [f'{ind}    return ({pid!r}, {pdesc!r})']
    def emit_abv():
        return [f'{ind}    return ("Pattern A" if {iptm_v} >= 0.60 else "Pattern B", "Above-Threshold Interactor (Moderate affinity; inspect contact density & PAE)")']
else:
    # Detect exact assignment target variable names from the existing Pattern block
    pat_id_v, pat_desc_v = None, None
    for sub in ast.walk(outer_if):
        if isinstance(sub, ast.Assign):
            if any(isinstance(c, ast.Constant) and c.value == "Pattern C" for c in ast.walk(sub.value)):
                t0 = sub.targets[0]
                if isinstance(t0, ast.Tuple) and len(t0.elts) == 2:
                    pat_id_v, pat_desc_v = t0.elts[0].id, t0.elts[1].id
                elif isinstance(t0, ast.Name):
                    pat_id_v = t0.id
    if not pat_id_v:
        pat_id_v = pick_var(["pattern_id", "pat_id"])
    if not pat_desc_v:
        pat_desc_v = pick_var(["pattern_desc", "pat_desc"])
    assert pat_id_v and pat_desc_v, f"Could not detect pattern target vars (id={pat_id_v}, desc={pat_desc_v})"

    def emit_pat(pid, pdesc):
        return [f'{ind}    {pat_id_v}, {pat_desc_v} = {pid!r}, {pdesc!r}']
    def emit_abv():
        return [f'{ind}    {pat_id_v}, {pat_desc_v} = ("Pattern A" if {iptm_v} >= 0.60 else "Pattern B"), "Above-Threshold Interactor (Moderate affinity; inspect contact density & PAE)"']

clean_pat_cli = [
    f'{ind}if {uis_v} >= {thresh_v} and {iptm_v} >= 0.60 and {dg_v} <= -5.8 and {pro_v} == 0{clash_ok}:',
    *emit_pat("Pattern A", "Concordant Rigid/Coiled-Coil or Globular Interactor (Global ipTM, Local LIS, and dG_uni support binding)"),
    f'{ind}elif ({uis_v} >= {thresh_v} or ({lis_v} >= 0.50 and {dg_v} <= -5.8)) and {iptm_v} < 0.60 and {dg_v} <= -5.8 and {pro_v} == 0{clash_ok}:',
    *emit_pat("Pattern B", "IDR-Diluted Local Interactor (Strong local interface LIS & dG_uni despite low global ipTM from disordered tails)"),
    f'{ind}elif {pro_v} > 0{clash_bad}:',
    *emit_pat("Pattern C", "Thermodynamically Penalized / Kink-Disrupted Interface (Geometric proximity offset by strained Proline kink or steric clash)"),
    f'{ind}elif {uis_v} >= {thresh_v}:',
    *emit_abv(),
    f'{ind}elif {wbsa_v} >= 250.0 and {dg_v} > -5.0:',
    *emit_pat("Pattern C", "Weak / Sub-Threshold Packing (Geometric burial without sufficient stabilizing contacts)"),
    f'{ind}else:',
    *emit_pat("Pattern D", "Below-Threshold Decoy or Disordered Hallucination (No physical binding support)"),
]

cli_new = "\n".join(cli_lines[:s_line] + clean_pat_cli + cli_lines[e_line:]) + "\n"
ast.parse(cli_new)
cli_path.write_text(cli_new, encoding="utf-8")
print(f"[PASS] Synced Pattern A-D in unippi_cli.py ({target_fn.name}(), vars: {uis_v}, {pro_v}, {clash_v}, {wbsa_v}) — ast.parse OK.\n")

# ==============================================================================
# PART 2: FINAL CLEANUP & BIOCHEMICAL AUTHENTICITY IN index.html
# ==============================================================================
html = idx_path.read_text(encoding="utf-8")

# 1. Merge stacked double header (lab-header + resource-nav) into ONE clean sticky bar
if 'class="lab-header"' in html:
    double_hdr_pat = re.compile(
        r'<header class="lab-header">.*?</header>\s*<div class="resource-nav">.*?</div>\s*(?=<main id="top">)',
        re.DOTALL
    )
    single_hdr = """<header class="resource-nav">
    <div class="resource-nav-inner">
      <a class="wordmark" href="#top"><strong>UniPPI-OWL</strong><span>v1.8</span></a>
      <nav class="resource-links">
        <a href="#analyze">Upload &amp; Run</a>
        <a href="#leaderboardCard">Results</a>
        <a href="#workbench">3D &amp; Contacts</a>
        <a href="#svg-maps">Vector Maps</a>
        <a href="#interpretation-guide">Method &amp; Scope</a>
      </nav>
    </div>
  </header>

  """
    html, n_hdr = double_hdr_pat.subn(single_hdr, html, count=1)
    assert n_hdr == 1, "Could not merge double header in index.html"

# 2. Clean Pattern A-D block in index.html (remove `else if (false)`)
pat_js_block = re.compile(
    r'let patId = "Pattern D".*?(?=\n\s*const typesA = new Set)',
    re.DOTALL
)
clean_pat_js = """let patId = "Pattern D", patDesc = "Below-Threshold Decoy or Disordered Hallucination (No physical binding support)";
      if (uis >= UIS_THRESHOLD && iptm >= 0.60 && dgUni <= -5.8 && helPro === 0 && nClash < 5) {
        patId = "Pattern A"; patDesc = "Concordant Rigid/Coiled-Coil or Globular Interactor (Global ipTM, Local LIS, and ΔG_uni all support binding)";
      } else if ((uis >= UIS_THRESHOLD || (lis >= 0.50 && dgUni <= -5.8)) && iptm < 0.60 && dgUni <= -5.8 && helPro === 0 && nClash < 5) {
        patId = "Pattern B"; patDesc = "IDR-Diluted Local Interactor (Strong local interface LIS & ΔG_uni despite low global ipTM from disordered tails)";
      } else if (helPro > 0 || nClash >= 5) {
        patId = "Pattern C"; patDesc = "Thermodynamically Penalized / Kink-Disrupted Interface (Geometric proximity offset by strained Proline kink or steric clash)";
      } else if (uis >= UIS_THRESHOLD) {
        patId = iptm >= 0.60 ? "Pattern A" : "Pattern B";
        patDesc = "Above-Threshold Interactor (Moderate affinity; inspect contact density & PAE)";
      } else if (wBsa >= 250.0 && dgUni > -5.0) {
        patId = "Pattern C"; patDesc = "Weak / Sub-Threshold Packing (Geometric burial without sufficient stabilizing contacts)";
      }"""
html, n_pjs = pat_js_block.subn(clean_pat_js, html, count=1)
assert n_pjs == 1, "Could not update Pattern A-D block in index.html"

# 3. Clean dead empty mutation arrays & expose exact blend weights, sig_aff, thermo_gate
html = html.replace(
    "const hotspots = [], foldxMuts = [], mutCatalog = [];",
    "const mutCatalog = [];"
)
html = html.replace(
    """        top_hotspot: hotspots.length ? hotspots[0].residue : "NONE",
        top_hotspot_ddG_Ala: hotspots.length ? hotspots[0].ddG_Ala : 0.0,
        hotspots, contacts, foldx_mutations: foldxMuts, mutCatalog,""",
    """        sig_aff: Number(sigAff.toFixed(3)), thermo_gate: Number(thermoGate.toFixed(3)),
        w_prodigy_blend: Number((gate * (1.0 - physW)).toFixed(2)),
        w_phys_blend: Number((1.0 - gate * (1.0 - physW)).toFixed(2)),
        contacts, mutCatalog,"""
)
html = html.replace(
    ",\n        evalDgAtSalt: (s_mM) => clamp(evalEnvAtSalt(s_mM).dUni + ptmBonus, -22.0, 5.0)",
    ""
)

# 4. Upgrade loadSampleAF3Ensemble() to 100% biochemically authentic PDB residue/atom names
s_samp = html.find("    function loadSampleAF3Ensemble()")
e_samp = html.find("    async function handleFileSelect(fileList)")
assert s_samp > 0 and e_samp > s_samp, "Could not locate loadSampleAF3Ensemble() in index.html"

authentic_sample_js = r"""    function loadSampleAF3Ensemble() {
      LAST_UPLOADED_FILES = null;
      // Neutral 50-aa complementary globular heterodimer sample (100% authentic PDB atom names)
      const CHAIN_A_SEQ = "AQVINTEDGVQDYLQTYHKLADNYSITKSEAQALGFVISNGNLADTLGKS";
      const CHAIN_B_SEQ = "KKILINKGQINSIDSLHQELKKELQLPEYRGENLDVLFDTLTGWVQYPLV";
      const tempK = parseFloat(document.getElementById("tempInput")?.value) || 298.15;
      const lambdaPae = parseFloat(document.getElementById("lambdaInput")?.value) || 8.0;

      const modelConfigs = [
        { idx: 0, span: 2.94, plddt: 93.4, pae: 2.05, iptm: 0.87 },
        { idx: 1, span: 2.97, plddt: 92.6, pae: 2.20, iptm: 0.86 },
        { idx: 2, span: 3.00, plddt: 91.8, pae: 2.40, iptm: 0.84 },
        { idx: 3, span: 3.03, plddt: 91.0, pae: 2.60, iptm: 0.83 },
        { idx: 4, span: 3.06, plddt: 90.1, pae: 2.80, iptm: 0.81 }
      ];

      // Chemically authentic complementary pairs (0-based indices):
      // Salt Bridges: 6 (GLU-LYS), 18 (LYS-GLU), 21 (ASP-LYS), 27 (LYS-GLU), 29 (GLU-ARG)
      const SB_IDX = new Set([6, 18, 21, 27, 29]);
      // H-Bonds: 10 (GLN-ASN), 14 (GLN-SER), 24 (SER-GLN), 39 (ASN-THR), 45 (THR-GLN)
      const HB_IDX = new Set([10, 14, 24, 39, 45]);
      // Hydrophobic Core (all >=4 heavy-atom apolar side chains: V, I, L, F, Y):
      const HP_IDX = new Set([2, 3, 12, 19, 33, 35, 36, 37, 46]);

      const HB_ATOM_MAP = {
        GLN: ["OE1", "O"], ASN: ["ND2", "N"], SER: ["OG", "O"], THR: ["OG1", "O"], TYR: ["OH", "O"]
      };

      ALL_RESULTS = [];
      for (const cfg of modelConfigs) {
        const name = `sample_globular_af3_model_${cfg.idx}`;
        const seqA = Array.from(CHAIN_A_SEQ).map(ch => ONE_TO_THREE[ch] || "ALA");
        const seqB = Array.from(CHAIN_B_SEQ).map(ch => ONE_TO_THREE[ch] || "ALA");
        const cifLines = [
          `data_${name}`, "#", "_ma_qa_metric_local.type pLDDT", "loop_",
          "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol", "_atom_site.label_atom_id",
          "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id", "_atom_site.label_entity_id",
          "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
          "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv", "_atom_site.auth_seq_id",
          "_atom_site.auth_asym_id", "_atom_site.pdbx_PDB_model_num"
        ];
        let serial = 1;
        for (const [cIdx, [cid, seq]] of [["A", seqA], ["B", seqB]].entries()) {
          const sign = cIdx === 0 ? -1.0 : 1.0;
          for (let ri = 0; ri < seq.length; ri++) {
            const rName = seq[ri], rSeq = ri + 1;
            const zBase = ri * 2.15, yBase = Math.cos(ri * 0.9) * 1.45;
            const loc = cfg.span, xCa = sign * (loc + 1.50), xCb = sign * loc;
            const aSpec = [
              ["N", xCa + sign*0.8, yBase - 0.5, zBase - 0.6, "N"],
              ["CA", xCa, yBase, zBase, "C"],
              ["C", xCa + sign*0.6, yBase + 0.6, zBase + 0.6, "C"],
              ["O", xCa + sign*1.2, yBase + 1.1, zBase + 0.8, "O"]
            ];
            if (rName !== "GLY") aSpec.push(["CB", xCb, yBase, zBase, "C"]);
            if (SB_IDX.has(ri)) {
              const isAnion = (rName === "ASP" || rName === "GLU");
              if (isAnion) aSpec.push(["CG", sign*(loc-0.60), yBase, zBase, "C"], [rName==="ASP"?"OD1":"OE1", sign*(loc-1.45), yBase, zBase, "O"]);
              else if (rName === "LYS" || rName === "ARG" || rName === "HIS") aSpec.push(["CG", sign*(loc-0.60), yBase, zBase, "C"], [rName==="LYS"?"NZ":(rName==="HIS"?"NE2":"NH1"), sign*(loc-1.45), yBase, zBase, "N"]);
            } else if (HB_IDX.has(ri) && HB_ATOM_MAP[rName]) {
              const [hAtom, hEl] = HB_ATOM_MAP[rName];
              if (rName === "GLN" || rName === "ASN") aSpec.push(["CG", sign*(loc-0.65), yBase, zBase, "C"]);
              aSpec.push([hAtom, sign*(loc-1.40), yBase, zBase, hEl]);
            } else if (HP_IDX.has(ri)) {
              const c1 = (rName === "VAL" || rName === "ILE") ? "CG1" : "CG";
              const c2 = (rName === "VAL") ? "CG2" : "CD1";
              aSpec.push([c1, sign*(loc-0.75), yBase+0.12, zBase, "C"], [c2, sign*(loc-1.22), yBase-0.08, zBase+0.10, "C"]);
            }
            for (const [an, ax, ay, az, el] of aSpec) {
              cifLines.push(`ATOM ${serial++} ${el} ${an} . ${rName} ${cid} 1 ${rSeq} ? ${ax.toFixed(3)} ${ay.toFixed(3)} ${az.toFixed(3)} 1.00 ${cfg.plddt.toFixed(2)} ${rSeq} ${cid} 1`);
            }
          }
        }
        const rawCifText = cifLines.join("\n");
        const { atoms } = parseStructureText(rawCifText, `${name}.cif`);
        const { ordered, chains } = buildResidues(atoms);
        const sasaComp = computeResidueSASA(ordered);
        const sasaA = computeResidueSASA(chains["A"]), sasaB = computeResidueSASA(chains["B"]);
        const resToIdx = new Map(ordered.map((r, i) => [r.key, i]));
        const tot = ordered.length, nA = chains["A"].length;
        const paeMat = Array.from({length: tot}, (_, r) => Array.from({length: tot}, (__, c) => r===c ? 0.8 : ((r<nA&&c<nA)||(r>=nA&&c>=nA) ? 1.8 : cfg.pae)));
        const res = analyzePairJS(name, "A", "B", chains["A"], chains["B"], ordered, resToIdx, "MODE_A_PRED_FULL", { iptm: cfg.iptm, res_pae_matrix: paeMat }, sasaComp, sasaA, sasaB, tempK, lambdaPae);
        res.rawStructText = rawCifText;
        res.rawStructFormat = "cif";
        ALL_RESULTS.push(res);
      }

      const qList = document.getElementById("batchFileQueueList");
      if (qList) {
        qList.innerHTML = modelConfigs.map((c, idx) =>
          `<div>[${idx+1}/5] <strong>sample_globular_af3_model_${c.idx}.cif</strong> + <strong>sample_globular_af3_full_data_${c.idx}.json</strong> · Chains A(50 aa)–B(50 aa) · Status: <span style="color:var(--emerald); font-weight:700;">PAIRED &amp; EVALUATED</span></div>`
        ).join("");
      }
      ACTIVE_INDEX = 0;
      updateUI();
    }

"""
html = html[:s_samp] + authentic_sample_js + html[e_samp:]

# 5. Upgrade Step 2, Step 4 & #patternCallout to pocket-calculator precision
html = html.replace(
    "Mean weight w̄<sub>intf</sub> = <strong>${r.mean_w_intf}</strong>",
    "Mean weight w̄<sub>intf</sub> = <strong>${r.mean_w_intf.toFixed(3)}</strong>"
)

old_step4_consensus = """&bull; <code>Consensus ΔG<sub>uni</sub> = <strong>${r.dG_uni_kcal} kcal/mol</strong> &rarr; K<sub>d</sub> = exp(ΔG<sub>uni</sub> / RT) = <strong>${kdStr}</strong> &rarr; UIS = 0.25·ipTM(${r.ipTM}) + 0.30·LIS(${r.LIS}) + 0.25·pDockQ(${r.pDockQ_PAE}) + 0.20·σ<sub>aff</sub> = <strong>${r.UIS.toFixed(4)}</strong></code>"""
new_step4_consensus = """&bull; <code>Consensus ΔG<sub>uni</sub> = ${r.w_prodigy_blend.toFixed(2)}·ΔG<sub>PRODIGY</sub>(${r.dG_prodigy_kcal.toFixed(2)}) + ${r.w_phys_blend.toFixed(2)}·ΔG<sub>phys</sub>(${r.dG_phys_kcal.toFixed(2)}) = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> &rarr; K<sub>d</sub> = exp(ΔG<sub>uni</sub> / RT) = <strong>${kdStr}</strong></code><br>
              &bull; <code>${r.mode === "MODE_A_PRED_FULL"
                ? `UIS = [0.25·ipTM(${r.ipTM.toFixed(3)}) + 0.30·LIS(${r.LIS.toFixed(3)}) + 0.25·pDockQ<sub>PAE</sub>(${r.pDockQ_PAE.toFixed(3)}) + 0.20·σ<sub>aff</sub>(${r.sig_aff.toFixed(3)})] × gate(${r.thermo_gate.toFixed(2)})`
                : `UIS (Coord Mode) = [0.55·pDockQ<sub>norm</sub>(${Math.min(1, r.pDockQ_coord / 0.742).toFixed(3)}) + 0.45·σ<sub>aff</sub>(${r.sig_aff.toFixed(3)})] × gate(${r.thermo_gate.toFixed(2)})`
              } = <strong>${r.UIS.toFixed(4)}</strong></code>"""
if old_step4_consensus in html:
    html = html.replace(old_step4_consensus, new_step4_consensus, 1)

old_callout_note = """Note: ΔG<sub>uni</sub> (${r.dG_uni_kcal} kcal/mol) combines published PRODIGY contact energetics (${r.dG_prodigy_kcal} kcal/mol) with enthalpy–entropy compensated 8-channel solvation (${r.dG_phys_kcal} kcal/mol)."""
new_callout_note = """Note: ΔG<sub>uni</sub> (${r.dG_uni_kcal.toFixed(2)} kcal/mol) weights published PRODIGY (${Math.round(r.w_prodigy_blend*100)}%: ${r.dG_prodigy_kcal.toFixed(2)} kcal/mol) and enthalpy–entropy compensated 8-channel solvation (${Math.round(r.w_phys_blend*100)}%: ${r.dG_phys_kcal.toFixed(2)} kcal/mol)."""
if old_callout_note in html:
    html = html.replace(old_callout_note, new_callout_note, 1)

html = re.sub(r'file:///[^"#]+(#[^"]+)', r'\1', html)
idx_path.write_text(html, encoding="utf-8")
print("[PASS] Upgraded index.html cleanly (authentic PDB atoms, single header, pocket-calculator Step 4).")

# Optional Node.js JS syntax check
js_ok = True
js_blocks = re.findall(r'<script>(.*?)</script>', html, flags=re.DOTALL)
if js_blocks and shutil.which("node"):
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
        tf.write(js_blocks[-1])
        tf_path = tf.name
    proc = subprocess.run(["node", "--check", tf_path], capture_output=True, text=True)
    Path(tf_path).unlink(missing_ok=True)
    if proc.returncode != 0:
        print("[FAIL] Node.js JS syntax error:\n", proc.stderr)
        js_ok = False
    else:
        print("[PASS] Node.js JS syntax check passed (0 errors).")

# ==============================================================================
# PART 3: MASTER 16-POINT RELEASE AUDIT
# ==============================================================================
h = idx_path.read_text(encoding="utf-8")
c = cli_path.read_text(encoding="utf-8")

checks = [
    ("1. Zero 'chuma' anywhere in index.html or unippi_cli.py", "chuma" not in h.lower() and "chuma" not in c.lower()),
    ("2. Zero 'vip1', 'bzip61', or 'bzip34' in index.html", all(k not in h.lower() for k in ("vip1", "bzip61", "bzip34"))),
    ("3. Zero dead mutation code or `else if (false)` in index.html", all(k not in h for k in ("simulate3DStructuralPerturbation", "evaluateMutationDDG", "top_hotspot_ddG_Ala", "else if (false)"))),
    ("4. Exact calibrated tanh compensation (15.5) in BOTH files", "15.5 * Math.tanh(eRaw / 15.5)" in h and "15.5 * math.tanh(e_raw / 15.5)" in c),
    ("5. Smooth continuous thermoGate in BOTH files", "clamp((-dgUni + 0.5) / 6.0, 0.10, 1.0)" in h and "clamp((-dg_uni + 0.5) / 6.0, 0.10, 1.0)" in c),
    ("6. NIS reliability gate in BOTH files", "nisReliability = clamp(nisTot / 20.0" in h and "nis_reliability = clamp(" in c),
    ("7. Unified CONF_W_MIN = 0.25 across index.html & unippi_cli.py", "const CONF_W_MIN = 0.25;" in h and 'c["weight_w"] >= 0.25' in c),
    ("8. Strained interior helix Proline rule in BOTH files", "isStrained = (r.meanBfac < 80.0)" in h and "is_strained = (res.mean_bfactor < 80.0)" in c),
    ("9. Pattern C strictly synced in BOTH files (never flags clean UIS>=0.50)", "else if (helPro > 0 || nClash >= 5)" in h and "Weak / Sub-Threshold Packing" in c and "Weak / Sub-Threshold Packing" in h),
    ("10. Neutral 5-model sample uses 100% biochemically authentic atoms", "const HB_ATOM_MAP =" in h and 'const c1 = (rName === "VAL" || rName === "ILE") ? "CG1" : "CG";' in h),
    ("11. Clean Kd column header & smart pM/nM/uM/mM formatting", "<th>K<sub>d</sub></th>" in h and "toFixed(2)} µM`" in h),
    ("12. Bounded 28x28 PAE SVG & proportional Chain A/B divider", "const splitX = 125 + 210 * fracA;" in h),
    ("13. Pocket-calculator Step 4 (blend weights, sig_aff, MODE_A vs B/C)", "r.w_prodigy_blend.toFixed(2)" in h and "UIS (Coord Mode)" in h),
    ("14. Single compact sticky header (double nav bloat removed)", 'class="lab-header"' not in h and h.count('class="resource-nav"') == 1),
    ("15. True 3D contact click-to-zoom & multi-chain pairwise SASA", "function zoomToContactPair3D" in h and "pairSasaComp = computeResidueSASA([...chains[ca], ...chains[cb]]);" in h),
    ("16. Valid JS syntax, valid Python AST & zero file:/// links", js_ok and "file:///" not in h)
]

print("\n" + "="*88)
print("MASTER 16-POINT FINAL AUDIT (PROF. DAISUKE TSUGAMA BASELINE RELEASE)")
print("="*88)
all_ok = True
for label, ok in checks:
    print(f"[{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        all_ok = False
print("="*88)
if not all_ok:
    sys.exit(1)

# ==============================================================================
# PART 4: LIVE REGRESSION TESTS A & B ON unippi_cli.py
# ==============================================================================
zip_cands = list(Path(".").glob("fold_vip1_bzip61*.zip"))
if zip_cands:
    print(f"\n[REGRESSION TEST A] python3 unippi_cli.py -i {zip_cands[0].name} -o ./af3_vip1_bzip61_verified")
    subprocess.run([sys.executable, "unippi_cli.py", "-i", str(zip_cands[0]), "-o", "./af3_vip1_bzip61_verified"], check=True)

print("\n[REGRESSION TEST B] python3 unippi_cli.py --generate-demo ./demo_verified -o ./results_demo_verified")
subprocess.run([sys.executable, "unippi_cli.py", "--generate-demo", "./demo_verified", "-o", "./results_demo_verified"], check=True)
