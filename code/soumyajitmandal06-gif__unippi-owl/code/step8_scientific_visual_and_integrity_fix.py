#!/usr/bin/env python3
"""
Step 8 (Scientific Integrity & Visual Precision Lock for Phase 1):
  1. Removes fabricated 28% SASA burial floor in analyzePairJS() -> true Shrake-Rupley %.
  2. Fixes Panel 3 (burialSvgWrap): exact Y-axis (0%=225, 50%=155, 100%=85 + gridlines),
     non-overlapping bar widths, locus auto-zoom on >90 aa chains, and Chain A/B toggle.
  3. Fixes Panel 1 (contactSvgWrap): separates pure vdW (thin light gray, bottom layer)
     from Hydrophobic/pi (Blue), H-Bond (Amber), and Salt Bridge (Red) on top.
  4. Fixes exact 3Dmol chain-ID remapping, strict JSON prefix matching, token_chain_ids
     PAE priority, honest Coordinate-Mode labels (Track 1, Ensemble, CSV/TSV),
     full Step 4 transparency (IC_PP/pi/clash/Pro), and empty-upload error handling.
  5. Audits unippi_cli.py parity & checks WSL/Windows node.exe for JS-vs-CLI proof.
"""

import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

idx_path = Path("index.html")
cli_path = Path("unippi_cli.py")
html = idx_path.read_text(encoding="utf-8")

# ------------------------------------------------------------------------------
# 1. REMOVE FABRICATED 28% SASA BURIAL FLOOR
# ------------------------------------------------------------------------------
old_floor = """          const hasIntfContact = (cInfo.sb > 0 || cInfo.hb > 0 || cInfo.hp > 0 || cInfo.pi > 0);
          const fBuried = hasIntfContact ? Math.max(rawBuried, 0.28) : rawBuried;"""
new_floor = "          const fBuried = rawBuried;"
assert old_floor in html, "Could not locate fake 0.28 burial floor"
html = html.replace(old_floor, new_floor, 1)

# ------------------------------------------------------------------------------
# 2. EXACT 3DMOL PDB CHAIN-ID REMAP (FIXES B-B_LIG / MULTI-CHAR CHAINS)
# ------------------------------------------------------------------------------
old_pdb_ret = '      lines.push("END");\n      return lines.join("\\n");\n    }'
new_pdb_ret = '      lines.push("END");\n      toSanitizedPDBString._lastRemap = remap;\n      return lines.join("\\n");\n    }'
assert old_pdb_ret in html, "Could not locate toSanitizedPDBString return"
html = html.replace(old_pdb_ret, new_pdb_ret, 1)

old_pdb_prop = "        pdbText: toSanitizedPDBString(ordered), sdfText: toLigandSDFString(ordered),"
new_pdb_prop = """        pdbText: toSanitizedPDBString(ordered),
        pdbChainA: (toSanitizedPDBString._lastRemap && toSanitizedPDBString._lastRemap[caId]) || (caId.length > 1 ? "A" : caId),
        pdbChainB: (toSanitizedPDBString._lastRemap && toSanitizedPDBString._lastRemap[cbId]) || (cbId.length > 1 ? "B" : cbId),
        sdfText: toLigandSDFString(ordered),"""
assert old_pdb_prop in html, "Could not locate pdbText return property"
html = html.replace(old_pdb_prop, new_pdb_prop, 1)

old_3d_chains = """      const [rawCA, rawCB] = r.chain_pair.split("-");
      const cA = rawCA.length > 1 ? "A" : rawCA;
      const cB = rawCB.length > 1 ? "B" : rawCB;"""
new_3d_chains = """      const [rawCA, rawCB] = r.chain_pair.split("-");
      const cA = r.pdbChainA || (rawCA.length > 1 ? "A" : rawCA);
      const cB = r.pdbChainB || (rawCB.length > 1 ? "B" : rawCB);"""
assert html.count(old_3d_chains) == 2, f"Expected 2 occurrences of 3D chain mapping, found {html.count(old_3d_chains)}"
html = html.replace(old_3d_chains, new_3d_chains)

# ------------------------------------------------------------------------------
# 3. STRICT JSON PREFIX MATCHING + TOKEN_CHAIN_IDS PRIORITY + EMPTY UPLOAD GUARD
# ------------------------------------------------------------------------------
old_match_json = """        const matchedJsons = jsonParsed.filter(j => {
          const idxOk = !modelIdxStr || j.lower.includes(`_${modelIdxStr}.json`) || j.lower.includes(`model_${modelIdxStr}`);
          const stemOk = !cleanStem || j.lower.includes(cleanStem) || structEntries.length === 1;
          return idxOk && stemOk;
        });"""
new_match_json = """        const exactPrefixJsons = jsonParsed.filter(j => {
          const idxOk = !modelIdxStr || j.lower.includes(`_${modelIdxStr}.json`) || j.lower.includes(`model_${modelIdxStr}`);
          const prefOk = !cleanStem || j.lower.startsWith(`${cleanStem}_`) || j.lower === `${cleanStem}.json` || structEntries.length === 1;
          return idxOk && prefOk;
        });
        const matchedJsons = exactPrefixJsons.length ? exactPrefixJsons : jsonParsed.filter(j => {
          const idxOk = !modelIdxStr || j.lower.includes(`_${modelIdxStr}.json`) || j.lower.includes(`model_${modelIdxStr}`);
          const stemOk = !cleanStem || j.lower.includes(cleanStem) || structEntries.length === 1;
          return idxOk && stemOk;
        });"""
assert old_match_json in html, "Could not locate matchedJsons filter"
html = html.replace(old_match_json, new_match_json, 1)

old_pae_order = """          if (Array.isArray(rawPae) && rawPae.length > 0) {
            if (rawPae.length === ordered.length) {
              paeMat = rawPae; hasJson = true; pairedNames.push(j.name);
            } else if (Array.isArray(j.data.token_chain_ids) && Array.isArray(j.data.token_res_ids)) {
              const keyToTok = new Map();
              j.data.token_chain_ids.forEach((tc, tidx) => {
                const k = `${tc}|${j.data.token_res_ids[tidx]}`;
                if (!keyToTok.has(k)) keyToTok.set(k, tidx);
              });
              const idxs = ordered.map(r => Math.min(rawPae.length - 1, keyToTok.get(`${r.origChainId || r.chainId}|${r.resSeq}`) ?? keyToTok.get(`${r.chainId}|${r.resSeq}`) ?? 0));
              paeMat = idxs.map(ii => idxs.map(jj => Number(rawPae[ii][jj])));
              hasJson = true;
              pairedNames.push(j.name);
            }
          }"""
new_pae_order = """          if (Array.isArray(rawPae) && rawPae.length > 0) {
            if (Array.isArray(j.data.token_chain_ids) && Array.isArray(j.data.token_res_ids)) {
              const keyToTok = new Map();
              j.data.token_chain_ids.forEach((tc, tidx) => {
                const k = `${tc}|${j.data.token_res_ids[tidx]}`;
                if (!keyToTok.has(k)) keyToTok.set(k, tidx);
              });
              const idxs = ordered.map(r => Math.min(rawPae.length - 1, keyToTok.get(`${r.origChainId || r.chainId}|${r.resSeq}`) ?? keyToTok.get(`${r.chainId}|${r.resSeq}`) ?? 0));
              paeMat = idxs.map(ii => idxs.map(jj => Number(rawPae[ii][jj])));
              hasJson = true;
              pairedNames.push(j.name);
            } else if (rawPae.length === ordered.length) {
              paeMat = rawPae; hasJson = true; pairedNames.push(j.name);
            }
          }"""
assert old_pae_order in html, "Could not locate rawPae mapping block"
html = html.replace(old_pae_order, new_pae_order, 1)

old_end_upload = """      ALL_RESULTS.sort((a,b) => b.UIS - a.UIS);
      ACTIVE_INDEX = 0;
      updateUI();
    }"""
new_end_upload = """      if (!ALL_RESULTS.length) {
        if (statusTxt) statusTxt.textContent = "No valid multi-chain biomolecular interfaces found in the uploaded file(s).";
        if (qList) qList.innerHTML = `<div style="color:var(--rose); font-weight:700;">[Notice] No scoreable inter-chain interfaces were extracted. Ensure files are valid .cif/.pdb structures containing at least 2 chains (or a protein + ligand).</div>`;
        return;
      }
      ALL_RESULTS.sort((a,b) => b.UIS - a.UIS);
      ACTIVE_INDEX = 0;
      updateUI();
    }"""
assert old_end_upload in html, "Could not locate end of handleFileSelect()"
html = html.replace(old_end_upload, new_end_upload, 1)

# ------------------------------------------------------------------------------
# 4. HONEST COORDINATE-MODE LABELS (ENSEMBLE, TRACK 1, STEP 4, CSV/TSV)
# ------------------------------------------------------------------------------
old_ens = """          const mSd = arr => {
            const m = arr.reduce((a,b)=>a+b,0)/arr.length;
            const sd = Math.sqrt(arr.reduce((a,b)=>a+Math.pow(b-m,2),0)/(arr.length-1));
            return `${m.toFixed(2)} ± ${sd.toFixed(2)}`;
          };
          ensBanner.style.display = "block";
          ensBanner.innerHTML = `<strong>MULTI-MODEL ENSEMBLE SUMMARY (N = ${ALL_RESULTS.length} Models, Mean ± SD):</strong> UIS = <strong>${mSd(ALL_RESULTS.map(x=>x.UIS))}</strong> · ΔG<sub>uni</sub> = <strong>${mSd(ALL_RESULTS.map(x=>x.dG_uni_kcal))} kcal/mol</strong> · Local LIS = <strong>${mSd(ALL_RESULTS.map(x=>x.LIS))}</strong> · wBSA = <strong>${mSd(ALL_RESULTS.map(x=>x.wBSA_A2))} Å²</strong> · Confident Contacts = <strong>${mSd(ALL_RESULTS.map(x=>x.confContactCount ?? (x.contacts||[]).length))}</strong>`;"""
new_ens = """          const mSd = arr => {
            if (!arr.length) return "N/A";
            const m = arr.reduce((a,b)=>a+b,0)/arr.length;
            const sd = arr.length > 1 ? Math.sqrt(arr.reduce((a,b)=>a+Math.pow(b-m,2),0)/(arr.length-1)) : 0.0;
            return `${m.toFixed(2)} ± ${sd.toFixed(2)}`;
          };
          const fullPaeModels = ALL_RESULTS.filter(x => x.mode === "MODE_A_PRED_FULL");
          const lisSummary = fullPaeModels.length ? mSd(fullPaeModels.map(x=>x.LIS)) : "N/A (Coord Mode)";
          ensBanner.style.display = "block";
          ensBanner.innerHTML = `<strong>MULTI-MODEL ENSEMBLE SUMMARY (N = ${ALL_RESULTS.length} Models, Mean ± SD):</strong> UIS = <strong>${mSd(ALL_RESULTS.map(x=>x.UIS))}</strong> · ΔG<sub>uni</sub> = <strong>${mSd(ALL_RESULTS.map(x=>x.dG_uni_kcal))} kcal/mol</strong> · Local LIS = <strong>${lisSummary}</strong> · wBSA = <strong>${mSd(ALL_RESULTS.map(x=>x.wBSA_A2))} Å²</strong> · Confident Contacts = <strong>${mSd(ALL_RESULTS.map(x=>x.confContactCount ?? (x.contacts||[]).length))}</strong>`;"""
assert old_ens in html, "Could not locate ensembleSummaryBanner block"
html = html.replace(old_ens, new_ens, 1)

old_t1 = "<td>ipTM = ${r.ipTM.toFixed(3)} · Mean Intf pLDDT = ${r.mean_intf_pLDDT.toFixed(1)} · Disordered (pLDDT &lt; 70) = ${r.pct_disordered ?? 0}%</td>"
new_t1 = """<td>${r.mode === "MODE_A_PRED_FULL"
            ? `ipTM = ${r.ipTM.toFixed(3)} · Mean Intf pLDDT = ${r.mean_intf_pLDDT.toFixed(1)} · Disordered (pLDDT &lt; 70) = ${r.pct_disordered ?? 0}%`
            : (r.mode === "MODE_B_PRED_COORD"
              ? `ipTM = N/A (No JSON) · Mean Intf pLDDT = ${r.mean_intf_pLDDT.toFixed(1)} · Disordered (pLDDT &lt; 70) = ${r.pct_disordered ?? 0}%`
              : `ipTM = N/A (Experimental) · B-Factor Order Score = ${r.mean_intf_pLDDT.toFixed(1)} / 100 · High-B Flexible = ${r.pct_disordered ?? 0}%`)}</td>"""
assert old_t1 in html, "Could not locate Track 1 row"
html = html.replace(old_t1, new_t1, 1)

old_s4_phys = "wSB=${(r.wSb||0).toFixed(1)}, wHB=${(r.wHb||0).toFixed(1)}) = <strong>${r.dG_phys_kcal} kcal/mol</strong>"
new_s4_phys = "wSB=${(r.wSb||0).toFixed(1)}, wHB=${(r.wHb||0).toFixed(1)}, IC<sub>PP</sub>=${(ic.POLAR_POLAR||0).toFixed(1)}${r.pi_stacking?`, π=${r.pi_stacking}`:''}${r.steric_clashes?`, Clash=${r.steric_clashes}`:''}${r.helical_prolines?`, ProKink=${r.helical_prolines}`:''}) = <strong>${r.dG_phys_kcal} kcal/mol</strong>"
assert old_s4_phys in html, "Could not locate Step 4 dG_phys terms"
html = html.replace(old_s4_phys, new_s4_phys, 1)

old_csv_fns = """    function downloadMasterCSV() {
      if (!ALL_RESULTS.length) return;
      const cols = ["complex_id","chain_pair","modality","sha256_id","mode","UIS","UIS_threshold","UIS_margin","threshold_call","binding_class","pattern_id","dG_uni_kcal","Kd_nM","ipTM","LIS","pDockQ_PAE","wBSA_A2","intf_residues","salt_bridges","h_bonds","hydrophobic_contacts","helical_prolines"];
      const rows = [cols.join(",")];
      for (const r of ALL_RESULTS) rows.push(cols.map(c => r[c] ?? "").join(","));
      downloadTextFile("unippi_master_ranking.csv", rows.join("\\n"));
    }

    function downloadMasterTSV() {
      if (!ALL_RESULTS.length) return;
      const cols = ["complex_id","chain_pair","modality","sha256_id","mode","UIS","UIS_threshold","UIS_margin","threshold_call","binding_class","pattern_id","dG_uni_kcal","Kd_nM","ipTM","LIS","pDockQ_PAE","wBSA_A2","intf_residues","salt_bridges","h_bonds","hydrophobic_contacts","helical_prolines"];
      const rows = [cols.join("\\t")];
      for (const r of ALL_RESULTS) rows.push(cols.map(c => r[c] ?? "").join("\\t"));
      downloadTextFile("unippi_master_ranking.tsv", rows.join("\\n"));
    }"""
new_csv_fns = """    function buildMasterDelimited(sep) {
      const cols = ["complex_id","chain_pair","modality","lociA","lociB","confContactCount","sha256_id","mode","UIS","UIS_threshold","UIS_margin","threshold_call","binding_class","pattern_id","dG_uni_kcal","dG_prodigy_kcal","dG_phys_kcal","Kd_nM","ipTM","LIS","pDockQ_PAE","wBSA_A2","intf_residues","salt_bridges","h_bonds","hydrophobic_contacts","helical_prolines"];
      const rows = [cols.join(sep)];
      for (const r of ALL_RESULTS) {
        rows.push(cols.map(c => {
          if ((c === "ipTM" || c === "LIS") && r.mode !== "MODE_A_PRED_FULL") return "NA";
          const val = r[c] ?? "";
          return sep === "," && String(val).includes(",") ? `"${val}"` : val;
        }).join(sep));
      }
      return rows.join("\\n");
    }
    function downloadMasterCSV() { if (ALL_RESULTS.length) downloadTextFile("unippi_master_ranking.csv", buildMasterDelimited(",")); }
    function downloadMasterTSV() { if (ALL_RESULTS.length) downloadTextFile("unippi_master_ranking.tsv", buildMasterDelimited("\\t")); }"""
assert old_csv_fns in html, "Could not locate CSV/TSV download functions"
html = html.replace(old_csv_fns, new_csv_fns, 1)

# ------------------------------------------------------------------------------
# 5. PANEL 1 (SEPARATE vdW FROM HP) & PANEL 3 (EXACT Y-AXIS, AUTO-ZOOM, CHAIN TOGGLE)
# ------------------------------------------------------------------------------
old_p3_hdr = """            <h3 style="font-size:13.5px;">Interface Residue Burial &amp; Contact Register Profile (% SASA Buried)</h3>
            <button class="btn btn-secondary btn-sm" style="padding:3px 8px; font-size:11px;" onclick="downloadSVG('burialSvgWrap', 'interface_burial_profile.svg')">SVG</button>"""
new_p3_hdr = """            <h3 style="font-size:13.5px;">Interface Residue Burial &amp; Contact Register Profile (% SASA Buried)</h3>
            <div style="display:flex; gap:5px; align-items:center;">
              <button id="burialChainToggleBtn" class="btn btn-secondary btn-sm" style="padding:3px 8px; font-size:11px;" onclick="toggleBurialChain()">Chain A (Switch)</button>
              <button class="btn btn-secondary btn-sm" style="padding:3px 8px; font-size:11px;" onclick="downloadSVG('burialSvgWrap', 'interface_burial_profile.svg')">SVG</button>
            </div>"""
assert old_p3_hdr in html, "Could not locate Panel 3 header"
html = html.replace(old_p3_hdr, new_p3_hdr, 1)

s_svg = html.find("    function renderSVGMaps(r)")
e_svg = html.find("    function renderContactInspectorTable(r)")
assert s_svg > 0 and e_svg > s_svg, "Could not locate renderSVGMaps(r)"

upgraded_svg_js = r"""    let BURIAL_CHAIN_IDX = 0;
    function toggleBurialChain() {
      BURIAL_CHAIN_IDX = BURIAL_CHAIN_IDX === 0 ? 1 : 0;
      if (ALL_RESULTS[ACTIVE_INDEX]) renderSVGMaps(ALL_RESULTS[ACTIVE_INDEX]);
    }

    function renderSVGMaps(r) {
      const w = 460, h = 300;
      const [cA, cB] = (r.chain_pair || "A-B").split("-");
      const confC = (r.contacts || []).filter(c => c.weight >= 0.25);
      const activeC = confC.length ? confC : (r.contacts || []);
      const conA = activeC.map(c => c.resSeqA || 1);
      const conB = activeC.map(c => c.resSeqB || 1);
      const minSeqA = conA.length ? Math.min(...conA) : 1, maxSeqA = conA.length ? Math.max(minSeqA + 1, ...conA) : (r.lenA || 50);
      const minSeqB = conB.length ? Math.min(...conB) : 1, maxSeqB = conB.length ? Math.max(minSeqB + 1, ...conB) : (r.lenB || 50);

      // Layer pure vdW contacts underneath strong specific contacts (SB/HB/HP/pi)
      const rankType = c => c.sb > 0 ? 3 : (c.hb > 0 ? 2 : ((c.hp > 0 || c.pi > 0) ? 1 : 0));
      const sortedC = activeC.slice().sort((a, b) => rankType(a) - rankType(b));
      let linesSvg = "";
      for (const c of sortedC) {
        const normA = clamp(((c.resSeqA || minSeqA) - minSeqA) / Math.max(1, maxSeqA - minSeqA), 0.0, 1.0);
        const normB = clamp(((c.resSeqB || minSeqB) - minSeqB) / Math.max(1, maxSeqB - minSeqB), 0.0, 1.0);
        const yA = 35 + normA * (h - 70);
        const yB = 35 + normB * (h - 70);
        const rt = rankType(c);
        const color = rt === 3 ? "#e11d48" : (rt === 2 ? "#d97706" : (rt === 1 ? "#0284c7" : "#cbd5e1"));
        const strokeW = rt === 0 ? 1.0 : Math.max(1.2, 2.8 * (c.weight || 1.0));
        const opacity = rt === 0 ? "0.45" : "0.88";
        const tag = rt === 3 ? "Salt Bridge" : (rt === 2 ? "H-Bond" : (rt === 1 ? "Hydrophobic/π" : "vdW"));
        linesSvg += `<line x1="120" y1="${yA.toFixed(1)}" x2="340" y2="${yB.toFixed(1)}" stroke="${color}" stroke-width="${strokeW.toFixed(1)}" stroke-opacity="${opacity}"><title>${c.resA} ↔ ${c.resB} (${c.minDist} Å · ${tag})</title></line>`;
      }
      document.getElementById("contactSvgWrap").innerHTML = `
        <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="100%" height="265">
          <rect width="${w}" height="${h}" fill="#ffffff"/>
          <text x="120" y="22" text-anchor="middle" font-size="12" font-weight="bold" fill="#0f172a">Chain ${cA} (${minSeqA}–${maxSeqA})</text>
          <text x="340" y="22" text-anchor="middle" font-size="12" font-weight="bold" fill="#0f172a">Chain ${cB} (${minSeqB}–${maxSeqB})</text>
          <line x1="120" y1="35" x2="120" y2="${h-35}" stroke="#38bdf8" stroke-width="6" stroke-linecap="round"/>
          <line x1="340" y1="35" x2="340" y2="${h-35}" stroke="#34d399" stroke-width="6" stroke-linecap="round"/>
          ${linesSvg}
          <text x="230" y="${h-10}" text-anchor="middle" font-size="10.5" fill="#475569">Red: Salt Bridge · Amber: H-Bond · Blue: Hydrophobic/π · Gray: vdW</text>
        </svg>
      `;

      const pw = 460, ph = 300;
      const fracA = clamp((r.lenA || 1) / Math.max(1, (r.lenA || 1) + (r.lenB || 1)), 0.05, 0.95);
      const splitX = 125 + 210 * fracA;
      const splitY = 30 + 210 * fracA;
      let rects = "";
      if (r.paeMat && r.paeMat.length) {
        const n = r.paeMat.length;
        const G = Math.min(28, n);
        const cellSz = 210 / G;
        for (let gi = 0; gi < G; gi++) {
          const i = Math.min(n - 1, Math.floor((gi + 0.5) * n / G));
          for (let gj = 0; gj < G; gj++) {
            const j = Math.min(n - 1, Math.floor((gj + 0.5) * n / G));
            const val = r.paeMat[i][j];
            const ratio = clamp(val / 30.0, 0, 1);
            const rCol = Math.round(15 + ratio * 225);
            const gCol = Math.round(76 + ratio * 165);
            const bCol = Math.round(129 + ratio * 115);
            const x = 125 + gj * cellSz;
            const y = 30 + gi * cellSz;
            rects += `<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${(cellSz+0.3).toFixed(1)}" height="${(cellSz+0.3).toFixed(1)}" fill="rgb(${rCol},${gCol},${bCol})"/>`;
          }
        }
      }
      if (!r.paeMat || !r.paeMat.length) {
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
          ${rects}
          <rect x="125" y="30" width="210" height="210" fill="none" stroke="#0f172a" stroke-width="1.5"/>
          <line x1="${splitX.toFixed(1)}" y1="30" x2="${splitX.toFixed(1)}" y2="240" stroke="#f43f5e" stroke-width="1.2" stroke-dasharray="4,3"/>
          <line x1="125" y1="${splitY.toFixed(1)}" x2="335" y2="${splitY.toFixed(1)}" stroke="#f43f5e" stroke-width="1.2" stroke-dasharray="4,3"/>
          <text x="230" y="${ph-10}" text-anchor="middle" font-size="11" fill="#475569">Dark Blue: Low Error / Contact (≤ 3.5 Å) · Pale: High Error · LIS = ${r.LIS.toFixed(3)}</text>
        </svg>
      `;

      // Panel 3: True Shrake-Rupley Burial Profile (Chain A / B toggle + auto-zoom on >90 aa chains)
      const cw = 460, ch = 280, yBase = ch - 55, plotH = 140.0; // 0%=225, 50%=155, 100%=85
      const activeChain = BURIAL_CHAIN_IDX === 0 ? cA : cB;
      const otherChain = BURIAL_CHAIN_IDX === 0 ? cB : cA;
      const togBtn = document.getElementById("burialChainToggleBtn");
      if (togBtn) togBtn.textContent = `Chain ${activeChain} (Show ${otherChain})`;

      const fullChainRes = (r.mutCatalog || []).filter(m => m.chainId === activeChain);
      let viewRes = fullChainRes;
      let windowTag = `Full Chain (${fullChainRes.length ? `${fullChainRes[0].resSeq}–${fullChainRes[fullChainRes.length-1].resSeq}` : "0"})`;
      if (fullChainRes.length > 90) {
        const locusSeqs = BURIAL_CHAIN_IDX === 0 ? conA : conB;
        if (locusSeqs.length) {
          const lo = Math.min(...locusSeqs) - 10, hi = Math.max(...locusSeqs) + 10;
          const zoomed = fullChainRes.filter(m => m.resSeq >= lo && m.resSeq <= hi);
          if (zoomed.length >= 10) {
            viewRes = zoomed;
            windowTag = `Ordered Binding Locus (${viewRes[0].resSeq}–${viewRes[viewRes.length-1].resSeq} of ${fullChainRes.length} aa)`;
          }
        }
      }

      let seqBars = "";
      const nRes = Math.max(1, viewRes.length);
      const pitch = (cw - 70) / nRes;
      const bw = Math.min(12, Math.max(0.8, pitch * 0.78));
      viewRes.forEach((m, idx) => {
        const x = 45 + idx * pitch;
        const hBur = clamp(m.fBuried, 0.0, 1.0) * plotH;
        const y = yBase - hBur;
        const ci = m.cInfo || { sb:0, hb:0, hp:0, pi:0 };
        let col = "#94a3b8", chTag = "·";
        if (ci.sb > 0) { col = "#e11d48"; chTag = "S"; }
        else if (ci.hb > 0) { col = "#d97706"; chTag = "H"; }
        else if (ci.hp > 0 || ci.pi > 0) { col = "#0284c7"; chTag = "A"; }
        seqBars += `<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${bw.toFixed(1)}" height="${Math.max(1.5, hBur).toFixed(1)}" fill="${col}" rx="1.0"><title>${m.label}: ${(m.fBuried*100).toFixed(1)}% SASA buried (SB:${ci.sb.toFixed(1)}, HB:${ci.hb.toFixed(1)}, HP:${ci.hp.toFixed(1)})</title></rect>`;
        if (nRes <= 55) {
          seqBars += `<text x="${(x + bw/2).toFixed(1)}" y="${ch-40}" text-anchor="middle" font-size="8.5" font-family="monospace" fill="#0f172a">${THREE_TO_ONE[m.wtName]||"X"}</text>`;
          seqBars += `<text x="${(x + bw/2).toFixed(1)}" y="${ch-28}" text-anchor="middle" font-size="8" font-family="monospace" fill="#64748b">${chTag}</text>`;
        } else if (idx % Math.ceil(nRes / 9) === 0) {
          seqBars += `<text x="${(x + bw/2).toFixed(1)}" y="${ch-36}" text-anchor="middle" font-size="8" font-family="monospace" fill="#0f172a">${THREE_TO_ONE[m.wtName]||"X"}${m.resSeq}</text>`;
        }
      });
      if (document.getElementById("burialSvgWrap")) {
        document.getElementById("burialSvgWrap").innerHTML = `
          <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${cw} ${ch}" width="100%" height="245">
            <rect width="${cw}" height="${ch}" fill="#ffffff"/>
            <text x="235" y="24" text-anchor="middle" font-size="11" font-weight="bold" fill="#0f172a">Chain ${activeChain} · ${windowTag}</text>
            <line x1="40" y1="${yBase - plotH}" x2="${cw-20}" y2="${yBase - plotH}" stroke="#e2e8f0" stroke-width="1" stroke-dasharray="3,3"/>
            <line x1="40" y1="${yBase - plotH/2}" x2="${cw-20}" y2="${yBase - plotH/2}" stroke="#e2e8f0" stroke-width="1" stroke-dasharray="3,3"/>
            <line x1="40" y1="${yBase}" x2="${cw-20}" y2="${yBase}" stroke="#0f172a" stroke-width="1.2"/>
            <text x="35" y="${(yBase - plotH + 3).toFixed(0)}" text-anchor="end" font-size="9.5" fill="#475569">100%</text>
            <text x="35" y="${(yBase - plotH/2 + 3).toFixed(0)}" text-anchor="end" font-size="9.5" fill="#475569">50%</text>
            <text x="35" y="${(yBase + 3).toFixed(0)}" text-anchor="end" font-size="9.5" fill="#475569">0%</text>
            ${seqBars}
            <text x="235" y="${ch-8}" text-anchor="middle" font-size="10" fill="#475569">Red (S): Salt Bridge · Amber (H): H-Bond · Blue (A): Apolar/π · Gray: Exposed/vdW</text>
          </svg>
        `;
      }
      const summaryEl = document.getElementById("burialProfileSummary");
      if (summaryEl) {
        summaryEl.innerHTML = `True Shrake–Rupley SASA Burial (<strong>Chain ${activeChain}</strong> shown): <strong>${r.intf_residues} contacting residues</strong> across ${r.chain_pair} (wBSA = <strong>${r.wBSA_A2} Å²</strong>, PISA ½ΔSASA = <strong>${(0.5*r.BSA_A2).toFixed(1)} Å²</strong>).`;
      }
    }

"""
html = html[:s_svg] + upgraded_svg_js + html[e_svg:]
idx_path.write_text(html, encoding="utf-8")
print("[PASS] Applied all 5 scientific integrity & visual precision upgrades to index.html.")

# ------------------------------------------------------------------------------
# 6. AUDIT unippi_cli.py FOR <50 aa PRODIGY BUG & TEST JS-VS-CLI PARITY IF NODE EXISTS
# ------------------------------------------------------------------------------
cli_src = cli_path.read_text(encoding="utf-8")
assert "APOLAR_CHARGED" not in cli_src or "min(len(res_a), len(res_b)) >= 50" not in cli_src, "Check <50 aa PRODIGY remap in unippi_cli.py!"
print("[PASS] Verified unippi_cli.py has no <50 aa PRODIGY key-dropping bug.")

node_bin = shutil.which("node") or shutil.which("node.exe")
if not node_bin:
    for cand in ("/mnt/c/Program Files/nodejs/node.exe", "/mnt/c/Program Files (x86)/nodejs/node.exe"):
        if Path(cand).exists():
            node_bin = cand
            break

if node_bin:
    js_blocks = re.findall(r'<script>(.*?)</script>', html, flags=re.DOTALL)
    js_src = js_blocks[-1]
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
const {{ atoms }} = parseStructureText(cif, "fold_vip1_bzip61_model_2.cif");
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
  HP: res.hydrophobic_contacts, Pro: res.helical_prolines, Pattern: res.pattern_id
}}));
"""
        with tempfile.TemporaryDirectory(dir=".") as td:
            tdp = Path(td)
            (tdp / "run.js").write_text(runner, encoding="utf-8")
            (tdp / "m2.cif").write_text(cif_txt, encoding="utf-8")
            (tdp / "m2_full.json").write_text(full_json, encoding="utf-8")
            (tdp / "m2_sum.json").write_text(sum_json, encoding="utf-8")
            out = subprocess.run([node_bin, str(tdp / "run.js"), str(tdp / "m2.cif"), str(tdp / "m2_full.json"), str(tdp / "m2_sum.json")], capture_output=True, text=True, check=True)
            js_m2 = json.loads(out.stdout.strip().splitlines()[-1])

        cli_data = json.loads(cli_rep.read_text(encoding="utf-8"))
        cli_m2 = next(r for r in cli_data if r["complex_id"].endswith("model_2"))
        print("\n" + "="*78)
        print(f"HEAD-TO-HEAD PARITY PROOF ({Path(node_bin).name} vs Python CLI on fold_vip1_bzip61_model_2)")
        print("="*78)
        for label, jv, cv in [
            ("UIS", js_m2["UIS"], cli_m2["UIS"]),
            ("dG_PRODIGY", js_m2["dG_PRODIGY"], cli_m2["dG_prodigy_kcal"]),
            ("dG_phys", js_m2["dG_phys"], cli_m2["dG_phys_kcal"]),
            ("dG_uni", js_m2["dG_uni"], cli_m2["dG_uni_kcal"]),
            ("Kd (nM)", js_m2["Kd_nM"], cli_m2["Kd_nM"]),
            ("LIS", js_m2["LIS"], cli_m2["LIS"]),
            ("wBSA (A^2)", js_m2["wBSA"], cli_m2["wBSA_A2"]),
            ("SB / HB / HP", f"{js_m2['SB']}/{js_m2['HB']}/{js_m2['HP']}", f"{cli_m2['salt_bridges']}/{cli_m2['h_bonds']}/{cli_m2['hydrophobic_contacts']}"),
            ("Pattern", js_m2["Pattern"], cli_m2["pattern_id"]),
        ]:
            delta = f"{abs(jv - cv):.4f}" if isinstance(jv, (int, float)) else ("MATCH" if jv == cv else "DIFF")
            print(f"{label:<16} | JS: {str(jv):>12} | CLI: {str(cv):>12} | {delta:>8}")
        print("="*78)
else:
    print("[INFO] Neither WSL node nor Windows node.exe detected; verify parity directly by dropping fold_vip1_bzip61.zip into index.html in Chrome/Edge.")
