"""Nearest BioEmu-1 training sequence for every SAXS construct (all 439), by component.

Search: MMseqs2 easy-search of the constructs against each component (sensitivity 7.5, E <= 10, up to 4,000
hits per query, full Smith-Waterman alignment of every hit, alignment-mode 3). Identity is MMseqs2 fident, the
identical positions over the alignment length; coverage is the share of the construct inside the alignment.
For each construct and component the nearest neighbour is the hit with the highest identity among those
covering at least half of the construct (also kept: the best identity at any coverage).

Components of BioEmu-1's training (Lewis et al. 2025): afdb (cluster input sequences) and afdb_member (all
structure representatives), cath1, cath2, megasim, megascale (curated wild types, a superset of the 361
used), other_md. desres is searched and reported separately (not in BioEmu-1's training). Octapeptides
(8 residues) cannot cover half of any construct and are not searched.
The nearest training identity of a construct is the maximum over afdb_member, cath1, cath2, megasim,
megascale and other_md.

Also per construct: whether its UniProtKB accession is an AFDB cluster member or input, and whether a
CATH superfamily inside the construct occurs in the training data (any component, or the MD components).
The metapredict disorder of each construct (script 14) is added for the comparisons that use it.
Writes tables/training_nearest.csv and tables/mmseqs_nearest_<component>.m8.
"""
import subprocess

import pandas as pd
from common import MMSEQS, SAXS_TABLE, TAB, THREADS

TRAIN = ["afdb_member", "cath1", "cath2", "megasim", "megascale", "other_md"]
SEARCH = ["afdb", *TRAIN, "desres"]
MD = ["cath1", "cath2", "megasim"]
MIN_COV = 0.5
TMP = TAB / "mmseqs_tmp"

saxs = pd.read_csv(SAXS_TABLE)
q = TAB / "fasta" / "saxs_constructs.fasta"
q.write_text("".join(f">{lab}\n{s}\n" for lab, s in zip(saxs.label, saxs.sequence, strict=True)))
cols = ["query", "target", "fident", "alnlen", "qstart", "qend", "qlen", "tstart", "tend", "tlen", "evalue", "bits",
        "qcov", "tcov"]

out = saxs[["label", "length"]].rename(columns={"length": "construct_length"}).set_index("label")
for comp in SEARCH:
    m8 = TAB / f"mmseqs_nearest_{comp}.m8"
    if not m8.exists():
        subprocess.run([MMSEQS, "easy-search", str(q), str(TAB / "fasta" / f"{comp}.fasta"), str(m8), str(TMP / comp),
                        "-s", "7.5", "-e", "10", "--max-seqs", "4000", "--alignment-mode", "3", "--threads", str(THREADS),
                        "-v", "1", "--format-output", ",".join(cols)], check=True)
    h = pd.read_csv(m8, sep="\t", header=None, names=cols)
    h["system"] = h.target.str.split("|", n=2).str[2]
    anyc = h.sort_values(["fident", "bits"], ascending=False).drop_duplicates("query").set_index("query")
    cov = h[h.qcov >= MIN_COV].sort_values(["fident", "bits"], ascending=False).drop_duplicates("query").set_index("query")
    out[f"id_{comp}"] = cov.fident.reindex(out.index).fillna(0.0)
    out[f"cov_{comp}"] = cov.qcov.reindex(out.index)
    out[f"hit_{comp}"] = cov.system.reindex(out.index)
    out[f"id_anycov_{comp}"] = anyc.fident.reindex(out.index).fillna(0.0)
    print(comp, "hits with coverage >= 0.5 for", int(cov.shape[0]), "constructs")

ids = out[[f"id_{c}" for c in TRAIN]]
out["nearest_identity"] = ids.max(axis=1)
out["nearest_component"] = ids.idxmax(axis=1).str.removeprefix("id_").where(out.nearest_identity > 0)
out["nearest_hit"] = [out.at[lab, f"hit_{c}"] if isinstance(c, str) else None
                      for lab, c in zip(out.index, out.nearest_component, strict=True)]
out["nearest_identity_md"] = out[[f"id_{c}" for c in MD]].max(axis=1)

# accession and superfamily overlap
mp = pd.read_csv(TAB / "per_protein_mapping.csv", keep_default_na=False, na_values=[""]).set_index("label")
mem = pd.read_csv(TAB / "afdb_members.csv.gz", usecols=["member", "input", "uniprot_entry"])
accs_member = set(mem.member) | set(mem.uniprot_entry.dropna())
accs_input = set(mem[mem.input].member) | set(mem[mem.input].uniprot_entry.dropna())
ta = pd.read_csv(TAB / "training_annotated.csv.gz", low_memory=False)
ta = ta[ta.component.isin(["afdb", *TRAIN])]


def sfset(s):
    return {x for x in str(s).split(";") if x and x != "nan"}


sf_any = set().union(*ta.cath_sf.map(sfset))
sf_md = set().union(*ta[ta.component.isin(MD)].cath_sf.map(sfset))
sf_afdb = set().union(*ta[ta.component == "afdb"].cath_sf.map(sfset))
out["uniprot"] = mp.uniprot.reindex(out.index)
out["clean"] = mp.clean.reindex(out.index).map({True: True, False: False, "True": True, "False": False})
out["disorder_class"] = mp.disorder_class.reindex(out.index)
out["uniprot_in_afdb_clusters"] = out.uniprot.map(lambda a: a in accs_member if isinstance(a, str) else None)
out["uniprot_is_afdb_input"] = out.uniprot.map(lambda a: a in accs_input if isinstance(a, str) else None)
sfs = mp.gene3d_superfamilies.reindex(out.index).map(sfset)
out["construct_superfamilies"] = sfs.map(lambda s: ";".join(sorted(s)))
out["n_superfamilies"] = sfs.map(len)
out["sf_in_training"] = [bool(s & sf_any) if s else None for s in sfs]
out["sf_in_md"] = [bool(s & sf_md) if s else None for s in sfs]
out["sf_in_afdb"] = [bool(s & sf_afdb) if s else None for s in sfs]
out["sf_group"] = ["none" if not s else ("present" if s & sf_any else "absent") for s in sfs]
out["sf_group_md"] = ["none" if not s else ("present" if s & sf_md else "absent") for s in sfs]
sm = pd.read_csv(TAB / "saxs_metapredict.csv").set_index("label")
out = out.join(sm[["dis_frac", "dis_longest", "dis_cat"]].add_prefix("mp_"))
out.reset_index().to_csv(TAB / "training_nearest.csv", index=False)
c = out[out.clean == True]  # noqa: E712
print(len(out), "constructs;", len(c), "clean")
for t in (0.3, 0.5, 0.9, 1.0):
    print(f"clean with nearest identity >= {t}: {(c.nearest_identity >= t).sum()}")
print(c.sf_group.value_counts().to_dict(), c.sf_group_md.value_counts().to_dict())
