"""Annotate every sequence set with the same three descriptors: CATH superfamilies, predicted disorder and
organism group.

CATH superfamilies
  SAXS constructs   CATH-Gene3D domains inside the construct (script 07)
  afdb, other_md    CATH-Gene3D cross-references of the UniProtKB entry (the whole chain is the training input)
  cath1, cath2      the CATH domain identifier of the simulated system, looked up in the CATH v4.4 domain list
                    (v4.3 list for identifiers retired since)
  megasim, megascale, desres
                    best MMseqs2 hit against CATH v4.4 S100 domain sequences (E <= 1e-3, at least 60% of the
                    query aligned); designed miniproteins mostly have none
  CATH reference    the S35 clusters of the domain list (script 08)
Disorder: metapredict 3.1 (Emenecker et al. 2021, Biophys J) on every sequence, the same predictor for every set: the share
  of residues scoring at least 0.5 and the longest run of such residues. Categories as in script 08: mostly
  disordered (>= 50% of residues), a long disordered region (>= 30 residues, under 50%), neither.
Organism group: UniProtKB lineage (afdb, other_md); source organism of the PDB entry (cath1, cath2 and the
  natural MEGAscale domains, named after a PDB entry); "Designed" for the de novo MEGAscale designs.
Writes tables/training_annotated.csv.gz, tables/saxs_metapredict.csv and tables/cath_s35_metapredict.csv.gz.
"""
import json
import os
import re
import subprocess

import importlib.metadata
import metapredict
import numpy as np
import pandas as pd
from common import DL, MMSEQS, SAXS_TABLE, TAB, THREADS, ORG_CATS, S, org_category, record_version

MMSEQS_RELEASE = os.environ.get("MMSEQS_RELEASE", "18-8cc5c")    # the release tag of the binary used
ROOT_TMP = TAB / "mmseqs_tmp"
ROOT_TMP.mkdir(exist_ok=True)
SUB = str.maketrans({"U": "C", "O": "K", "B": "N", "Z": "Q", "J": "L", "X": "G"})


def disorder(seqs: dict) -> pd.DataFrame:
    """metapredict per-residue scores -> share disordered, longest disordered run, category."""
    clean = {k: v.upper().translate(SUB) for k, v in seqs.items()}
    out = []
    keys = list(clean)
    for i in range(0, len(keys), 5000):
        res = metapredict.predict_disorder({k: clean[k] for k in keys[i:i + 5000]}, return_numpy=True)
        for k, (_, sc) in res.items():
            d = np.asarray(sc) >= 0.5
            run = best = 0
            for v in d:
                run = run + 1 if v else 0
                best = max(best, run)
            out.append((k, float(d.mean()), int(best)))
    df = pd.DataFrame(out, columns=["key", "dis_frac", "dis_longest"]).set_index("key")
    df["dis_cat"] = np.where(df.dis_frac >= 0.5, "mostly disordered",
                             np.where(df.dis_longest >= 30, "long disordered region", "neither"))
    return df


def mmseqs_best_cath(fasta, name):
    """Best hit of each query against CATH v4.4 S100 domain sequences."""
    out = TAB / f"mmseqs_cath_{name}.m8"
    if not out.exists():
        subprocess.run([MMSEQS, "easy-search", str(fasta), str(DL / "cath" / "cath-domain-seqs-S100.fa"), str(out),
                        str(ROOT_TMP / name), "-s", "7.5", "-e", "1e-3", "--threads", str(THREADS), "-v", "1",
                        "--format-output", "query,target,fident,alnlen,qstart,qend,qlen,tlen,evalue,bits,qcov,tcov"],
                       check=True)
    m = pd.read_csv(out, sep="\t", header=None, names=["query", "target", "fident", "alnlen", "qstart", "qend", "qlen",
                                                       "tlen", "evalue", "bits", "qcov", "tcov"])
    m = m[m.qcov >= 0.6].sort_values("bits", ascending=False).drop_duplicates("query")
    m["domain"] = m.target.str.split("|").str[2].str.split("/").str[0]
    return m


# ------------------------------------------------------------------ CATH domain list (v4.4, v4.3 fallback)
def domain_sf(path):
    d = pd.read_csv(path, sep=r"\s+", comment="#", header=None, usecols=range(5), names=["domain", "C", "A", "T", "H"])
    return dict(zip(d.domain, d.C.astype(str) + "." + d.A.astype(str) + "." + d["T"].astype(str) + "." + d.H.astype(str),
                    strict=True))


sf44 = domain_sf(DL / "cath" / "cath-domain-list.txt")
p43 = DL / "cath" / "cath-domain-list-v4_3_0.txt"
if not p43.exists():
    p43.write_bytes(S.get("https://download.cathdb.info/cath/releases/all-releases/v4_3_0/cath-classification-data/"
                          "cath-domain-list-v4_3_0.txt", timeout=600).content)
sf43 = domain_sf(p43)

ts = pd.read_csv(TAB / "training_sequences.csv.gz", low_memory=False)
ts = ts[ts.component != "afdb_member"].reset_index(drop=True)
ts["key"] = ts.component + "|" + ts.index.astype(str)

# superfamilies
ts["cath_sf"] = ""
ts["cath_source"] = ""
g3 = ts.gene3d.fillna("").map(lambda s: ";".join(sorted({x.strip() for x in s.split(";") if x.strip()})))
k = ts.component.isin(["afdb"])
ts.loc[k, "cath_sf"], ts.loc[k, "cath_source"] = g3[k], "Gene3D (UniProtKB)"
for i in ts.index[ts.component == "other_md"]:
    acc = ts.at[i, "uniprot"]
    e = json.loads((DL / "bioemu_training" / "other_md" / f"{acc}.json").read_text())
    ids = sorted({x["id"].replace("G3DSA:", "") for x in e.get("uniProtKBCrossReferences", []) if x["database"] == "Gene3D"})
    ts.at[i, "cath_sf"], ts.at[i, "cath_source"] = ";".join(ids), "Gene3D (UniProtKB)"
for i in ts.index[ts.component.isin(["cath1", "cath2"])]:
    dom = ts.at[i, "cath_domain"]
    ts.at[i, "cath_sf"] = sf44.get(dom) or sf43.get(dom, "")
    ts.at[i, "cath_source"] = "CATH v4.4 domain list" if dom in sf44 else ("CATH v4.3 domain list" if dom in sf43 else "not found")
q = ts[ts.component.isin(["megasim", "megascale", "desres"])]
fa = TAB / "fasta" / "small_domains_query.fasta"
fa.write_text("".join(f">{k_}\n{s_}\n" for k_, s_ in zip(q.key, q.sequence, strict=True)))
hits = mmseqs_best_cath(fa, "small_domains").set_index("query")
for i in q.index:
    key = ts.at[i, "key"]
    if key in hits.index:
        dom = hits.at[key, "domain"]
        ts.at[i, "cath_sf"] = sf44.get(dom, "")
        ts.at[i, "cath_source"] = f"MMseqs2 best hit {dom} ({hits.at[key, 'fident']:.2f} identity)"
    else:
        ts.at[i, "cath_source"] = "no CATH hit"

# disorder: training sets, SAXS constructs, CATH S35 reference
dis = disorder(dict(zip(ts.key, ts.sequence, strict=True)))
ts = ts.join(dis, on="key")
saxs = pd.read_csv(SAXS_TABLE)
sd = disorder(dict(zip(saxs.label, saxs.sequence, strict=True)))
sd.rename_axis("label").to_csv(TAB / "saxs_metapredict.csv")
cath_seqs = {}
for block in (DL / "cath" / "cath-domain-seqs-S35.fa").read_text().split(">")[1:]:
    head, *seq = block.splitlines()
    cath_seqs[head.split("|")[2].split("/")[0]] = "".join(seq)
cd = disorder(cath_seqs)
cd["length"] = [len(cath_seqs[k_]) for k_ in cd.index]
cd["sf"] = [sf44.get(k_, "") for k_ in cd.index]
cd.rename_axis("domain").to_csv(TAB / "cath_s35_metapredict.csv.gz", compression={"method": "gzip", "mtime": 0})
record_version("metapredict", {"version": importlib.metadata.version("metapredict"),
                               "source": "https://github.com/idptools/metapredict", "threshold": 0.5})
record_version("MMseqs2", {"version": f"release {MMSEQS_RELEASE}",
                           "commit": subprocess.run([MMSEQS, "version"], capture_output=True, text=True).stdout.strip(),
                           "source": "https://github.com/soedinglab/MMseqs2/releases (static linux-avx2 binary)"})

# organism groups
def lineage_names(s):
    return [re.sub(r"\s*\([^)]*\)$", "", x.strip()) for x in s.split(",")] if isinstance(s, str) else []


ts["org_group"] = None
k = ts.component.isin(["afdb", "other_md"])
ts.loc[k, "org_group"] = [org_category(t, n, lineage_names(lin)) for t, n, lin in
                          zip(ts.loc[k, "taxid"], ts.loc[k, "organism"], ts.loc[k, "lineage"], strict=True)]
pdb_of = {}
for i in ts.index[ts.component.isin(["cath1", "cath2", "megasim", "megascale", "desres"])]:
    s = ts.at[i, "system"]
    m = re.match(r"^(\d[0-9A-Za-z]{3})([A-Za-z0-9]?)", s)
    if ts.at[i, "component"] in ("cath1", "cath2"):
        pdb_of[i] = (s[:4].upper(), s[4])
    elif ts.at[i, "component"] == "desres":
        pdb_of[i] = (s.upper(), None)
    elif m and (len(s) == 4 or s[4] in "_."):
        pdb_of[i] = (s[:4].upper(), None)
    else:
        ts.at[i, "org_group"] = "Designed"
Q = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id polymer_entities{
  rcsb_polymer_entity_container_identifiers{ auth_asym_ids }
  rcsb_entity_source_organism{ ncbi_taxonomy_id ncbi_scientific_name taxonomy_lineage{ name } } } } }"""
cache = DL / "rcsb" / "training_entries.json"
ids = sorted({v[0] for v in pdb_of.values()})
got = json.loads(cache.read_text()) if cache.exists() else {}
todo = [x for x in ids if x not in got]
for j in range(0, len(todo), 200):
    r = S.post("https://data.rcsb.org/graphql", json={"query": Q, "variables": {"ids": todo[j:j + 200]}}, timeout=300)
    r.raise_for_status()
    for e in r.json()["data"]["entries"] or []:
        got[e["rcsb_id"]] = e
    got.update({x: None for x in todo[j:j + 200] if x not in got})
cache.write_text(json.dumps(got))
for i, (pdb, chain) in pdb_of.items():
    e = got.get(pdb)
    if not e:
        ts.at[i, "org_group"] = "Other"
        continue
    ents = e["polymer_entities"] or []
    ent = next((x for x in ents if chain and chain in (x["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] or [])),
               ents[0] if ents else None)
    src = (ent or {}).get("rcsb_entity_source_organism") or []
    if src:
        o = src[0]
        ts.at[i, "taxid"], ts.at[i, "organism"] = o.get("ncbi_taxonomy_id"), o.get("ncbi_scientific_name")
        ts.at[i, "org_group"] = org_category(o.get("ncbi_taxonomy_id"), o.get("ncbi_scientific_name") or "",
                                             [x["name"] for x in o.get("taxonomy_lineage") or []])
    else:
        ts.at[i, "org_group"] = "Other"
assert set(ts.org_group.dropna()) <= set(ORG_CATS) | {"Designed"}
ts.to_csv(TAB / "training_annotated.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})
print(ts.groupby("component").agg(n=("key", "size"), with_sf=("cath_sf", lambda s: int((s != "").sum())),
                                  mostly_dis=("dis_cat", lambda s: int((s == "mostly disordered").sum())),
                                  median_len=("length", "median")))
print(ts.groupby("component").org_group.value_counts().unstack(fill_value=0))
print("CATH MD domains not found:", ts[ts.cath_source == "not found"].system.tolist())
