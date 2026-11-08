"""extract_helpers.py: the few functions the paper's phosphorylation analyses use from the unpublished temperature project (BENDER-MT).

The analyses import salt_ph_phospho.py, which pulls in four more BENDER-MT files. Only the definitions they reach are released, each
copied verbatim from its source, so the rest of that project stays private. Starting from WANTED, the tool follows names inside each
definition (same-module names, and attributes of the module aliases) until nothing new turns up. Standard library only.
    python3 -I tools/extract_helpers.py   ->   helpers/{salt_ph_phospho,md_loader,fdt_params,basin_frames}.py, audit/helpers_extracted.txt
"""
import ast
from pathlib import Path

MT = Path("/Volumes/User Homes/Nu Project/BENDER-MT/code")
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "helpers"
# source module -> released module
MODS = {"salt_ph_phospho": "salt_ph_phospho", "TODO_T2_trajectory_loader": "md_loader", "TODO_T2_fdt_computation": "fdt_params",
        "basin_susceptibility": "basin_frames"}
# module aliases used inside the sources
ALIAS = {"L": "TODO_T2_trajectory_loader", "fdt": "TODO_T2_fdt_computation", "bs": "basin_susceptibility"}
WANTED = {"salt_ph_phospho": {"CAL", "BETA_C", "MAX_ELEMS", "site_matrix", "calvados_terms", "T", "bc"},
          "TODO_T2_trajectory_loader": {"K_FRAMES", "load_residue_params"},
          "TODO_T2_fdt_computation": {"EPS_AH", "RC_AH", "RC_YU", "calvados_residue_params", "get_dielectric_params", "pair_setup", "R_GAS"},
          "basin_susceptibility": {"frames_calvados"}}
STD = {"os", "sys", "numpy", "pandas", "scipy", "traceback", "argparse", "multiprocessing"}


def top_level(tree):
    """{name: node} for functions, classes and simple assignments at module level, including those inside a top-level if or try."""
    out = {}

    def add(body):
        for n in body:
            if isinstance(n, (ast.FunctionDef, ast.ClassDef)):
                out[n.name] = n
            elif isinstance(n, ast.Assign):
                for t in n.targets:
                    for x in ast.walk(t):
                        if isinstance(x, ast.Name):
                            out[x.id] = n
            elif isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name):
                out[n.target.id] = n
    add(tree.body)
    return out


def main():
    src = {m: (MT / f"{m}.py").read_text() for m in MODS}
    trees = {m: ast.parse(s) for m, s in src.items()}
    tops = {m: top_level(t) for m, t in trees.items()}
    keep = {m: {} for m in MODS}
    other_refs, imports = {}, {m: set() for m in MODS}
    for m, t in trees.items():
        for n in t.body:
            if isinstance(n, ast.Import):
                for a in n.names:
                    imports[m].add((a.name, a.asname))
            elif isinstance(n, ast.ImportFrom):
                for a in n.names:
                    imports[m].add((f"{n.module}.{a.name}", a.asname))
    todo = [(m, w) for m, ws in WANTED.items() for w in ws]
    used_mods = {m: set() for m in MODS}
    while todo:
        m, name = todo.pop()
        node = tops[m].get(name)
        assert node is not None, (m, name)
        if id(node) in {id(x) for x in keep[m].values()}:
            keep[m][name] = node
            continue
        keep[m][name] = node
        for x in ast.walk(node):
            if isinstance(x, ast.Attribute) and isinstance(x.value, ast.Name):
                base = x.value.id
                if base in ALIAS:
                    used_mods[m].add(base)
                    if x.attr != "bc":           # bc is bender_common, released in full
                        todo.append((ALIAS[base], x.attr))
                elif base in ("mh",) and base not in tops[m]:
                    other_refs.setdefault(f"{m}: {base}.{x.attr}", name)
            if isinstance(x, ast.Name) and isinstance(x.ctx, ast.Load) and x.id != name:
                if x.id in tops[m]:
                    todo.append((m, x.id))
                elif x.id in ALIAS:
                    used_mods[m].add(x.id)
    OUT.mkdir(exist_ok=True)
    log = []
    for m, rel in MODS.items():
        # source order, each statement once, text exactly as in the source
        # a module-level "if" that resets a kept name (a path fallback, say) belongs with it
        guards = [n for n in trees[m].body if isinstance(n, ast.If)
                  and {x.id for x in ast.walk(n) if isinstance(x, ast.Name) and isinstance(x.ctx, ast.Store)} & set(keep[m])]
        nodes = sorted({id(n): n for n in list(keep[m].values()) + guards}.values(), key=lambda n: n.lineno)
        names = [getattr(n, "name", None) or ast.get_source_segment(src[m], n).split("=")[0].strip() for n in nodes]
        head = [f'"""{rel}.py: definitions copied verbatim from the unpublished {m}.py (BENDER-MT) by tools/extract_helpers.py.\nDo not edit by hand."""']
        head += ["import os", "import sys", "", "import numpy as np", "import pandas as pd", "from scipy import sparse", ""]
        head += [ast.get_source_segment(src[m], n) for n in trees[m].body if isinstance(n, ast.ImportFrom) and n.module == "typing"]
        if m == "TODO_T2_trajectory_loader":
            head += ["sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', os.path.dirname(os.path.abspath(__file__))))",
                     "import bender_common as bc  # noqa: E402", ""]
        for alias in sorted(used_mods[m] | ({"L", "fdt", "bs"} if m == "salt_ph_phospho" else set())):
            if MODS[ALIAS[alias]] != rel:
                head.append(f"import {MODS[ALIAS[alias]]} as {alias}  # noqa: E402")
        lines = src[m].splitlines()
        body = []
        for n in nodes:
            start = min([n.lineno] + [d.lineno for d in getattr(n, "decorator_list", [])])
            # keep a comment block directly above a definition
            while start > 1 and lines[start - 2].startswith("#"):
                start -= 1
            body += lines[start - 1:n.end_lineno] + [""]
        (OUT / f"{rel}.py").write_text("\n".join(head + [""] + body))
        log.append(f"{rel}.py <- {m}.py: {len(nodes)} of {len(tops[m])} definitions: {', '.join(names)}")
    extra = [f"unresolved reference, {k}, used by {v}" for k, v in sorted(other_refs.items())]
    txt = "\n\n".join(log + extra) + "\n"
    (ROOT / "audit").mkdir(exist_ok=True)
    (ROOT / "audit/helpers_extracted.txt").write_text(txt)
    print(txt)


if __name__ == "__main__":
    main()
