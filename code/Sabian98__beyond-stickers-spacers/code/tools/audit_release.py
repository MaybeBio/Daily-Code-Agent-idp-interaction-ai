"""audit_release.py: what the code release of "Beyond Stickers and Spacers" must contain, and what is missing.

Starts from the figures the paper includes and the analyses it reports (SEEDS), then follows Python imports and references to
other scripts until nothing new turns up. Read-only on every source folder. Standard library only.
    python3 -I tools/audit_release.py   ->   audit/{manifest.csv, data_files.csv, open_questions.md}

manifest.csv     every script in the release: where it comes from, why it is included, its sha256, whether that hash matches
                 the last entry of gaivi/FROZEN.sha256 (the file as it was run), and how many absolute paths it holds
data_files.csv   every data file those scripts name, whether it exists here, and its size
open_questions   what a reader could not reproduce or what should not be published, for a person to settle
"""
import ast
import csv
import glob
import hashlib
import os
import re
from pathlib import Path

HOME = Path.home()
PAPER = HOME / "Downloads/Scripts/RECOMB_paper_template"
DYN = Path("/Volumes/User Homes/Nu Project/BENDER-DYN")
EV = DYN / "redo_evidence"
FIGS = DYN / "methods_paper/figures"
REAN = HOME / "Downloads/Scripts/bender_reanalysis"
OUT = Path(__file__).resolve().parent.parent / "audit"
SEARCH = [FIGS, EV, EV / "hub_predictor", EV / "gaivi", EV / "phospho_nmr_amend2", EV / "phospho_nmr_window", EV / "phospho_rauh",
          EV / "phospho_nmr_ceiling", REAN, REAN / "calvados", REAN / "mpipi", REAN / "bio_paper", Path("/Volumes/User Homes/Nu Project/BENDER-MT/code")]
# analyses the paper reports (Methods, Results, Supplementary Sections S1 to S17); work the paper does not use is left out on purpose.
SEEDS = """sep2_vs_sep3.py hub_frames.py gaivi/hub_labels_gaivi.py gaivi/submit_hub_labels.sh gaivi/verify_frozen.py
gaivi/design_sequences.py gaivi/design_glycine.py gaivi/analyze_glycine.py gaivi/design_arrangement_crossff.py
arrangement_test.py arrangement_test_amend_h4.py arrangement_test_degree.py fresh_ped_test.py ped_blind11_report.py
ped_extra_test.py pre_direct_check.py design_vs_experiment.py clinvar_hub.py clinvar_redo_predict.py clinvar_redo_test.py
llps_profiles.py llps_unevenness.py llps_unevenness_diagnostic.py llps_sticker_hub.py llps_mutations.py llps_csat*.py
llps_csat_degree.py llps_unevenness_degree.py llps_mutations_degree.py llps_degree_family.py degree_wetlab.py
hub_vs_degree_sim.py hub_vs_degree_models.py hub_vs_degree_iso.py hub_vs_degree_ped.py hub_vs_degree_ext.py hub_vs_degree_e49.py
hub_predictor/*.py phospho_network.py phospho_network_test.py phospho_network_posthoc.py phospho_direct_inputs.py
phospho_direct_inputs_amend1.py phospho_direct_labels.py phospho_direct_test.py phospho_direct_recharged_test.py
phospho_recharged_ess.py phospho_charge_test.py phospho_charge_posthoc_scd.py phospho_charge_posthoc_kappa.py
phospho_function_sites.py phospho_function_test.py phospho_nmr_inputs.py phospho_nmr_page4_prediction.py phospho_nmr_test.py
phospho_nmr_recharged_test.py phospho_nmr_rw_diagnostic.py phospho_nmr_amend2/make_inputs_nmr2.py
phospho_nmr_amend2/fetch_shifts_nmr2.py phospho_nmr_amend2/phospho_nmr_amend2_test.py phospho_nmr_window/phospho_nmr_window.py
phospho_nmr_ceiling/phospho_nmr_ceiling.py phospho_rauh/phospho_rauh.py phospho_rauh/make_inputs_rauh.py
phospho_rauh/phospho_rauh_direct.py recharged_last16/compare.py
agreement_ceiling.py arrangement_ff_table.py sim_null.py sim_null_v3.py sim_null_gg_pairs.py null_multi.py pitfall_theory.py
frames_reliability.py positional_null_check.py recompute_v3.sh label_recharged_v3.sh label_shuffles_v3.sh label_yf_v3.sh
update_recharged_v3_all.sh llps_degree_queue.sh ped_v2_build.py ped_v2_labels.py ped_v2_score.py ped_control.py
ped_residue_classes.py ped_restraint_information.py pre_reweight.py pre_forward_check.py hub_atlas.py consensus_hubs.py
context_stickers.py crossff_aa_degree.py structure_vs_agreement.py helicity_allatom.py dssp_numpy.py hub_overlap_by_ff.py
hub1000_analysis.py analyze_design.py clinvar_hub_conservation.py fetch_ped.sh fetch_ped_capped.sh fetch_ped_capped_multi.sh
hub_frames_summary.py hub_frames_crossff.py xlsx_to_csv.py""".split()
# steps the import trace cannot see: scripts that prepare what a figure script then reads
FIG_SEEDS = """case_fus_prep.py case_fus_raster.py figS_nmr_maps_prep.py figS_nmr_maps_prep2.py figS_reweighting_prep.py
fig2_crossff_export.py fig4_design_degree_export.py""".split()
SECRET = re.compile(r"(api[_-]?key|secret|passw(or)?d|token\s*=|hf_[A-Za-z0-9]{20,}|ghp_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9]{20,})", re.I)
EMAIL = re.compile(r"[\w.+-]+@[\w-]+\.[\w.]+")
ABS = re.compile(r"""['"](/(?:Users|Volumes|data|home|private)[^'"\n]*)['"]""")
DATA = re.compile(r"""['"]([^'"\n]*\.(?:csv|tsv|txt|pkl|npz|npy|json|fasta|xlsx|parquet|joblib))['"]""")
STDLIB_HINT = re.compile(r"^(numpy|pandas|scipy|sklearn|matplotlib|seaborn|mdtraj|networkx|joblib|openmm|pymupdf|fitz|PIL|requests|tqdm|lightgbm|xgboost|cugraph|finches|Bio|yaml|calvados|numba|statsmodels)$")


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def frozen():
    last = {}
    for line in open(EV / "gaivi/FROZEN.sha256", errors="ignore"):
        parts = line.split("#")[0].split()
        if len(parts) == 2:
            last[os.path.normpath(EV / "gaivi" / parts[1])] = parts[0]
    return last


def imports(path):
    try:
        tree = ast.parse(path.read_text(errors="ignore"))
    except SyntaxError:
        return None
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
            out.add(node.module.split(".")[0])
    return out


def resolve(name, near):
    for root in [near] + SEARCH:
        p = root / f"{name}.py"
        if p.exists():
            return p
    return None


def figure_scripts():
    """Each figure the paper includes, and the script that writes it (the script that names it in a savefig or output path)."""
    src = "".join(re.sub(r"(?<!\\)%.*", "", (PAPER / f).read_text()) for f in ("intro.tex", "results.tex", "appendix.tex"))
    names = sorted({Path(n).stem for n in re.findall(r"includegraphics(?:\[[^\]]*\])?\{figs/([^}]+)\}", src)})
    scripts = [Path(p) for p in glob.glob(str(FIGS / "*.py"))]
    texts = {p: p.read_text(errors="ignore") for p in scripts if not p.name.startswith("._")}
    found, missing = {}, []
    for n in names:
        hit = [p for p, t in texts.items() if n in t or p.stem == n]
        exact = [p for p in hit if p.stem == n]
        if exact or hit:
            found[n] = (exact or hit)[:2]
        else:
            missing.append(n)
    return names, found, missing


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frz = frozen()
    names, found, fig_missing = figure_scripts()
    todo = []
    for fig, ps in found.items():
        todo += [(p, f"figure {fig}") for p in ps]
    seed_missing = []
    for s in SEEDS:
        hits = [Path(p) for p in glob.glob(str(EV / s)) if not Path(p).name.startswith("._")]
        if not hits:
            seed_missing.append(s)
        todo += [(p, "analysis reported in the paper") for p in hits]
    todo += [(FIGS / s, "prepares a figure") for s in FIG_SEEDS if (FIGS / s).exists()]
    seen, unresolved, syntax = {}, {}, []
    data, external = {}, set()
    while todo:
        p, why = todo.pop()
        if p in seen:
            seen[p].add(why)
            continue
        seen[p] = {why}
        if p.suffix != ".py":
            continue
        mods = imports(p)
        if mods is None:
            syntax.append(p)
            continue
        text = p.read_text(errors="ignore")
        for m in mods:
            q = resolve(m, p.parent)
            if q:
                todo.append((q, f"imported by {p.name}"))
            elif STDLIB_HINT.match(m):
                external.add(m)
            elif m not in __import__("sys").stdlib_module_names:
                unresolved.setdefault(m, set()).add(p.name)
        for ref in re.findall(r"""['"/]([\w.-]+\.(?:py|sh))['"]""", text):
            q = resolve(ref[:-3], p.parent) if ref.endswith(".py") else (p.parent / ref if (p.parent / ref).exists() else None)
            if q and q != p:
                todo.append((q, f"called by {p.name}"))
        for ref in set(DATA.findall(text)):
            data.setdefault(ref, set()).add(p)
    rows, secrets, abspaths, oldhome = [], [], 0, []
    for p, why in sorted(seen.items(), key=lambda kv: str(kv[0])):
        text = p.read_text(errors="ignore")
        h = sha(p)
        fz = frz.get(os.path.normpath(p))
        state = "not in ledger" if fz is None else "matches ledger" if fz == h else "DIFFERS from ledger"
        n_abs = len(ABS.findall(text))
        abspaths += n_abs
        if "taseefrahman" in text:
            oldhome.append(p)
        for ln, line in enumerate(text.splitlines(), 1):
            if SECRET.search(line) or EMAIL.search(line):
                secrets.append((p, ln, line.strip()[:110]))
        rows.append({"file": str(p), "why": "; ".join(sorted(why))[:160], "sha256": h, "ledger": state, "absolute_paths": n_abs,
                     "bytes": p.stat().st_size})
    with open(OUT / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    drows = []
    for ref, users in sorted(data.items()):
        hits = [Path(ref)] if ref.startswith("/") else [u.parent / ref for u in users] + [EV / ref, FIGS / ref]
        ok = next((h for h in hits if h.exists()), None)
        drows.append({"named_as": ref, "exists": bool(ok), "bytes": ok.stat().st_size if ok else "", "resolved": str(ok or ""),
                      "used_by": " ".join(sorted(u.name for u in users))[:120]})
    with open(OUT / "data_files.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(drows[0]))
        w.writeheader()
        w.writerows(drows)
    differs = [r for r in rows if r["ledger"].startswith("DIFFERS")]
    notin = [r for r in rows if r["ledger"] == "not in ledger"]
    missing_data = [d for d in drows if not d["exists"]]
    size = sum(d["bytes"] for d in drows if d["exists"])
    L = ["# Open questions for the code release", "",
         f"Generated by tools/audit_release.py. Scripts in the release: {len(rows)}; figures in the paper: {len(names)}; data files named: {len(drows)} "
         f"({size / 1e9:.2f} GB found); absolute paths to rewrite: {abspaths}.", "",
         f"Scripts by ledger state: {sum(r['ledger'] == 'matches ledger' for r in rows)} match the hash ledger, {len(differs)} differ, {len(notin)} are not in it.",
         "", f"External packages used: {', '.join(sorted(external))}.", "",
         f"## Figures with no script found ({len(fig_missing)})", ""] + [f"- {n}" for n in fig_missing] + \
        ["", f"## Seeds that matched no file ({len(seed_missing)})", ""] + [f"- {s}" for s in seed_missing] + \
        ["", f"## Imports that resolve to nothing ({len(unresolved)})", ""] + [f"- `{m}`, in {', '.join(sorted(v))[:120]}" for m, v in sorted(unresolved.items())] + \
        ["", f"## Files that no longer parse ({len(syntax)})", ""] + [f"- {p}" for p in syntax] + \
        ["", f"## Scripts that differ from the hash ledger ({len(differs)})", ""] + [f"- {r['file']}" for r in differs] + \
        ["", f"## Possible secrets or personal addresses, check before publishing ({len(secrets)})", ""] + \
        [f"- {p.name}:{ln}: `{t}`" for p, ln, t in secrets[:60]] + \
        ["", f"## Scripts with the old home folder /Users/taseefrahman ({len(oldhome)})", ""] + [f"- {p}" for p in oldhome] + \
        ["", f"## Data files named by scripts but not found ({len(missing_data)}, first 60)", ""] + \
        [f"- `{d['named_as']}`, used by {d['used_by']}" for d in missing_data[:60]]
    (OUT / "open_questions.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[:9]))
    print(f"figures without script: {len(fig_missing)}, seeds missing: {len(seed_missing)}, unresolved imports: {len(unresolved)}, no parse: {len(syntax)}, "
          f"differ: {len(differs)}, secret hits: {len(secrets)}, old home: {len(oldhome)}, data missing: {len(missing_data)}")


if __name__ == "__main__":
    main()
