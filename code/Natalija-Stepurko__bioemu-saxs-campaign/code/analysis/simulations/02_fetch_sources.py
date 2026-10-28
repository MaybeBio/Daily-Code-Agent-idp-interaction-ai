"""Downloads at pinned versions: the CALVADOS 3 residue parameters, PULCHRA 3.06 (compiled here with gcc),
Pepsi-SAXS 3.0 (Linux binary, academic use; not redistributed), and the protein sets the CALVADOS parameters
were trained and validated on (design §7). Writes results/simulations/calvados_sets.csv (set, file, name,
sequence) and results/simulations/sources.json (URLs, commits, sha256 of every download)."""

import hashlib
import json
import subprocess
import tarfile
import zipfile

import pandas as pd
import requests
from common import (
    CALVADOS3_RESIDUES,
    CAO_COMMIT,
    DL,
    OUT,
    PAPERS_COMMIT,
    PEPSI,
    PEPSI_URL,
    PULCHRA,
    PULCHRA_URL,
    TOOLS,
)

PAPERS = f"https://raw.githubusercontent.com/KULL-Centre/papers/{PAPERS_COMMIT}"
CAO = f"https://raw.githubusercontent.com/KULL-Centre/_2024_Cao_CALVADOSCOM/{CAO_COMMIT}"
# (set, role, url): optimisation inputs of each published parameter set
SETS = [
    (
        "CALVADOS 1",
        "training",
        f"{PAPERS}/2021/CG-IDPs-Tesei-et-al/optimization/5_M1/o1_init_proteinsRgs.pkl",
    ),
    ("CALVADOS 1", "training", f"{PAPERS}/2021/CG-IDPs-Tesei-et-al/optimization/5_M1/o1_init_proteins.pkl"),
    (
        "CALVADOS 2",
        "training",
        f"{PAPERS}/2022/CG-cutoffs-Tesei-et-al/optimization/M1_20/1_init_proteinsRgs.pkl",
    ),
    (
        "CALVADOS 2",
        "training",
        f"{PAPERS}/2022/CG-cutoffs-Tesei-et-al/optimization/M1_20/1_init_proteins.pkl",
    ),
    ("CALVADOS 3", "training", f"{CAO}/data/IDPs_MDPsCOM_2.2_0.08_2/allproteins.pkl"),
    ("CALVADOS 3", "validation", f"{CAO}/data/IDPs_MDPsCOM_2.2_0.08_2_validate/allproteins_validate.pkl"),
    ("CALVADOS 3", "validation", f"{CAO}/data/IDPs_MDPsCOM_2.2_0.08_2_validate/IDPsRgs_validate.pkl"),
    ("CALVADOS 3", "validation", f"{CAO}/data/IDPs_MDPsCOM_2.2_0.08_2_validate/MultiDomainsRgs_validate.pkl"),
]


def fetch(url, path):
    if not path.exists():
        r = requests.get(url, timeout=300)
        r.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(r.content)
    return hashlib.sha256(path.read_bytes()).hexdigest()


src = {
    "calvados3_residues": {
        "url": CALVADOS3_RESIDUES,
        "sha256": fetch(CALVADOS3_RESIDUES, DL / "residues_CALVADOS3.csv"),
    }
}
# PULCHRA: MIT licence, C source
src["pulchra"] = {"url": PULCHRA_URL, "sha256": fetch(PULCHRA_URL, DL / "pulchra_306.tgz")}
if not PULCHRA.exists():
    with tarfile.open(DL / "pulchra_306.tgz") as t:
        t.extractall(TOOLS, filter="data")
    subprocess.run(
        ["gcc", "-O3", "-o", "pulchra", "pulchra.c", "pulchra_data.c", "-lm"],
        cwd=PULCHRA.parent,
        check=True,
        capture_output=True,
    )
# Pepsi-SAXS: binary for academic use, downloaded and kept local
src["pepsi_saxs"] = {"url": PEPSI_URL, "sha256": fetch(PEPSI_URL, DL / "Pepsi-SAXS-Linux.zip")}
if not PEPSI.exists():
    zipfile.ZipFile(DL / "Pepsi-SAXS-Linux.zip").extractall(PEPSI.parent)
    PEPSI.chmod(0o755)
ver = subprocess.run([str(PEPSI), "--version"], capture_output=True, text=True).stdout
src["pepsi_saxs"]["version"] = next(
    (ln.split("version:", 1)[1].strip() for ln in ver.splitlines() if "version:" in ln), None
)
rows = []
src["calvados_sets"] = []
for name, role, url in SETS:
    p = DL / "calvados_sets" / url.split("/", 6)[-1].replace("/", "__")
    sha = fetch(url, p)
    df = pd.read_pickle(p)
    src["calvados_sets"].append({"set": name, "role": role, "url": url, "sha256": sha, "n": int(len(df))})
    for prot, r in df.iterrows():
        rows.append(
            {
                "set": name,
                "role": role,
                "file": url.rsplit("/", 1)[-1],
                "name": prot,
                "sequence": "".join(r["fasta"]) if not isinstance(r["fasta"], str) else r["fasta"],
            }
        )
t = pd.DataFrame(rows)
t.to_csv(OUT / "calvados_sets.csv", index=False)
(OUT / "sources.json").write_text(json.dumps(src, indent=1) + "\n")
print(t.groupby(["set", "role"]).name.nunique().to_string(), "\ndistinct sequences:", t.sequence.nunique())
