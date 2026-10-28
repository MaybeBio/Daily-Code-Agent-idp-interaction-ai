"""One table of BioEmu-1 training sequences, by component, from the downloads of scripts 11 and 12.

Components (Lewis et al. 2025, Table S1 and Materials and Methods):
  afdb          AFDB pre-training clusters: the input sequence of each cluster (highest mean pLDDT member)
  afdb_member   every structure representative of those clusters (at least 80% identical to the input)
  cath1         ONE_cath1 MD, 50 CATH domains (sequence from the simulated topology)
  cath2         MSR_cath2 MD, CATH v4.3 domains of 50-200 residues (sequence from the simulated topology)
  megasim       MSR_megasim_merge MD, MEGAscale wild types (sequence from the simulated topology)
  megascale     MEGAscale curated wild types (Tsuboyama et al. 2023 Dataset2/Dataset3, mut_type "wt");
                the 361 used for fine-tuning are a subset that the paper does not list
  other_md      the other public MD sets of Table S1, as named proteins (full UniProtKB chains)
  desres        DESRES fast folders: not part of BioEmu-1's training (used for a separate model)
Counted without sequences: octapeptides (1100 systems of 8 residues) and MegaSim point mutants.
Writes tables/training_sequences.csv.gz and one FASTA per component in tables/fasta/.
"""
import json
import re

import pandas as pd
from common import BT, TAB

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H", "CYX": "C",
         "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K", "NLE": "L", "MSE": "M"}
Z = BT / "zenodo"
FA = TAB / "fasta"
FA.mkdir(exist_ok=True)


def topology_sequence(path):
    """One-letter sequence from the CA atoms of the first model of a topology.pdb; caps have no CA."""
    seq, seen = [], set()
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            key = (line[21], line[22:27])
            if key not in seen:
                seen.add(key)
                seq.append(THREE.get(line[17:20].strip(), "X"))
    return "".join(seq)


def count_systems(listing):
    """Number of system directories in a Zenodo archive listing (<dataset>/<system>/...)."""
    names = pd.read_csv(Z / "listings" / listing).name.str.split("/")
    return int(names[names.str.len() > 2].str[1].nunique())


rows = []
# AFDB
mem = pd.read_csv(TAB / "afdb_members.csv.gz")
for _, r in mem.iterrows():
    if not isinstance(r.sequence, str):
        continue
    base = {"system": r.member, "cluster": r.cluster, "sequence": r.sequence, "taxid": r.taxid,
            "organism": r.organism, "lineage": r.lineage, "gene3d": r.gene3d,
            "plddt_mean": r.plddt_mean}
    rows.append({"component": "afdb_member", **base})
    if r.input:
        rows.append({"component": "afdb", **base})
# MD from topologies
for comp, ds in (("cath1", "ONE_cath1"), ("cath2", "MSR_cath2"), ("megasim", "MSR_megasim_merge")):
    for d in sorted(p for p in (Z / "systems" / ds).iterdir() if p.is_dir()):
        name = d.name.split("_", 1)[1] if comp in ("cath1", "cath2") else d.name
        rows.append({"component": comp, "system": name, "sequence": topology_sequence(d / "topology.pdb"),
                     "cath_domain": name if comp in ("cath1", "cath2") else None})
# MEGAscale curated wild types
ms = pd.read_csv(BT / "megascale" / "Tsuboyama2023_Dataset2_Dataset3_20230416.csv.gz", low_memory=False)
wt = ms[ms.mut_type == "wt"].drop_duplicates("WT_name")
for _, r in wt.iterrows():
    rows.append({"component": "megascale", "system": r.WT_name.removesuffix(".pdb").replace("|", ":"),
                 "sequence": r.aa_seq})
# other public MD sets: named proteins, chains of polyproteins
named = json.loads((BT / "other_md" / "named_proteins.json").read_text())
for ds, prot, acc, chains in named:
    e = json.loads((BT / "other_md" / f"{acc}.json").read_text())
    full = e["sequence"]["value"]
    org = e["organism"]
    if chains:
        feats = [f for f in e.get("features", []) if f["type"] == "Chain"]
        for c in chains:      # "nsp3" ... "nsp16": UniProt names the chain "... nspN" or "Non-structural protein N"
            f = next(f for f in feats if re.search(rf"\b{c}$", f.get("description", ""))
                     or f.get("description") == f"Non-structural protein {c[3:]}")
            a, b = f["location"]["start"]["value"], f["location"]["end"]["value"]
            rows.append({"component": "other_md", "system": f"{ds}:{prot.split()[0]} {c}", "sequence": full[a - 1:b],
                         "uniprot": acc, "taxid": org["taxonId"], "organism": org["scientificName"],
                         "lineage": ", ".join(org.get("lineage", []))})
    else:
        rows.append({"component": "other_md", "system": f"{ds}:{prot}", "sequence": full, "uniprot": acc,
                     "taxid": org["taxonId"], "organism": org["scientificName"], "lineage": ", ".join(org.get("lineage", []))})
# DESRES fast folders
for p in sorted((BT / "desres").glob("*.fasta")):
    recs = [(r.splitlines()[0], "".join(r.splitlines()[1:])) for r in p.read_text().split(">") if r.strip()]
    head, seq = next((h, q) for h, q in recs if not set(q) <= set("ACGTUN"))     # the protein entity
    rows.append({"component": "desres", "system": p.stem, "sequence": seq, "pdb_title": head})

df = pd.DataFrame(rows)
df = df[df.sequence.str.len() > 0]
df["length"] = df.sequence.str.len()
# other_md: one row per distinct sequence (the spike appears in two Folding@home sets)
df = pd.concat([df[df.component != "other_md"], df[df.component == "other_md"].drop_duplicates("sequence")])
df.to_csv(TAB / "training_sequences.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
for comp, g in df.groupby("component"):
    with open(FA / f"{comp}.fasta", "w") as fh:
        for i, (s, q) in enumerate(zip(g.system, g.sequence, strict=True)):
            fh.write(f">{comp}|{i}|{s}\n{q}\n")
counts = {"octapeptides_systems": count_systems("ONE_octapeptides.zip.csv.gz"),
          "megasim_mutant_systems": count_systems("MSR_megasim_mutants_disp_allatom.zip.csv.gz"),
          "megascale_curated_rows": int(len(ms)), "megascale_curated_wild_types": int(len(wt)),
          "megascale_curated_unique_sequences": int(ms.aa_seq.nunique())}
(TAB / "training_counts.json").write_text(json.dumps(counts, indent=1))
print(df.groupby("component").length.describe().round(1))
print(counts)
