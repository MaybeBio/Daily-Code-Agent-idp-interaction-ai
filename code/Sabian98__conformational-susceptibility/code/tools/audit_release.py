"""audit_release.py: chooses the scripts of this release and writes audit/manifest.csv (every script of the project's code folder,
whether it is included, why, and its SHA-256). Included are the scripts behind the paper's figures, tables, and numbers, and every
script they import. Left out are slide-deck and infographic scripts, exploratory analyses the paper does not report, and internal
planning scripts. Standard library only.
    python3 -I tools/audit_release.py
"""
import ast
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CODE = Path("/Volumes/User Homes/Nu Project/BENDER-MT/code")
HELP = Path.home() / "Downloads/Scripts/bender_reanalysis"
EXTRA = ["gaivi_templates/run_calvados_cal.py", "gaivi_templates/run_openmpipi_cal.py", "gaivi_templates/submit_page4.sh",
         "submit_mpipi_recharged_fdt.sh", "submit_mpipi_temperature.sh"]
LEFT_OUT = {
    "slide deck or infographic, not in the paper": "deck_edit_solvation deck_edit_extensions deck_edit_mutation fig_infographic "
    "fig_arial_deckfigs fig_deck_v2 fig_total_response fig_composition_vs_compaction fig_leftover_breakdown fig_noise_floor "
    "fig_charge_compact_breakdown",
    "exploratory, not reported in the paper": "arg_fret_predictions bayes_fit bayes_proteome bremer_predict joint_fit own_hydrophobic "
    "tea_design total_response_basins pnt_salt_check residual_compactness_check page4_gamma_scan hub_temperature_pilot "
    "make_era_inputs era_cal_analysis",
    "file housekeeping for the simulation runs, not analysis": "ingest_temperature_reruns flag_truncated_temperatures",
    "internal planning for experiments, not in the paper": "fig_page4_nmr_brief page4_nmr_current page4_nmr_direct page4_effect_scan "
    "page4_mutation_scan page4_salt_scan make_page4_nmimic_inputs make_page4_5rq_inputs",
}
EXCL = {m: why for why, ms in LEFT_OUT.items() for m in ms.split()}


def imports(path, known):
    out = set()
    for n in ast.walk(ast.parse(path.read_text(errors="ignore"))):
        if isinstance(n, ast.Import):
            out |= {a.name.split(".")[0] for a in n.names}
        elif isinstance(n, ast.ImportFrom) and n.module:
            out.add(n.module.split(".")[0])
    return out & known


def readme_notes():
    notes = {}
    for line in (CODE / "README.md").read_text(errors="ignore").splitlines():
        if line.startswith("| `"):
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            for n in re.findall(r"`([A-Za-z0-9_]+)\.py`", cells[0]):
                notes[n] = re.sub(r"\s+", " ", cells[1])[:200] if len(cells) > 1 else ""
    return notes


def main():
    local = {p.stem: p for p in CODE.glob("*.py") if not p.name.startswith("._")}
    helpers = {p.stem: p for p in HELP.glob("*.py") if not p.name.startswith("._")}
    known = set(local) | set(helpers)
    inc = {m: "paper script" for m in local if m not in EXCL}
    stack = list(inc)
    while stack:                                   # every import of an included script is included as well
        m = stack.pop()
        p = local.get(m) or helpers[m]
        for d in imports(p, known):
            if d not in inc:
                inc[d] = f"imported by {m}"
                stack.append(d)
    notes = readme_notes()
    rows = []
    for m, p in sorted(local.items()):
        keep = m in inc
        rows.append([str(p), "yes" if keep else "no", inc[m] if keep else EXCL.get(m, "not used by the paper"), notes.get(m, ""),
                     hashlib.sha256(p.read_bytes()).hexdigest()])
    for m in sorted(set(inc) & set(helpers) - set(local)):
        p = helpers[m]
        rows.append([str(p), "yes", inc[m], "shared helper of the BENDER analysis code", hashlib.sha256(p.read_bytes()).hexdigest()])
    for e in EXTRA:
        p = CODE / e
        rows.append([str(p), "yes", "simulation template or job script", "", hashlib.sha256(p.read_bytes()).hexdigest()])
    with open(ROOT / "audit/manifest.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["file", "included", "reason", "note", "sha256"])
        w.writerows(rows)
    n = sum(r[1] == "yes" for r in rows)
    print(f"{n} files included, {len(rows) - n} left out")


if __name__ == "__main__":
    main()
