"""Composition of BioEmu-1's training data against the SAXS set and the CATH reference, and the overlap of
the SAXS constructs with the training data. Writes tables/training_summary.json and the comparison tables.

Sets
  SAXS        the 399 clean entries (constructs); CATH units are (UniProt accession, Gene3D superfamily)
              pairs inside the construct, as in script 08
  AFDB        the AFDB pre-training clusters, one input sequence per cluster; units are (cluster, Gene3D
              superfamily of the UniProtKB entry) pairs
  CATH MD     CATH1 and CATH2 simulated domains; one unit per domain (its CATH superfamily)
  MEGAscale   curated MEGAscale wild types and the MegaSim wild types, one row per distinct sequence;
              units from the best CATH hit (script 14)
  CATH S35    CATH v4.4 S35 clusters (composition, length and disorder of the S35 representative sequences)
  PDB 30%     RCSB 30% clusters of chains of at least 30 residues (organism only, as in script 08)
Disorder and organism are read from scripts 14 (training sets, SAXS constructs, CATH S35) and 07 (SAXS).
The MEGAscale domains called mostly disordered are split by origin (natural, de novo designed; Fisher exact
test) and their lengths compared with the rest (Mann-Whitney), in `megascale_mostly_disordered`.
Single proportions carry Wilson 95% intervals.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats
from common import DL, ORG_CATS, TAB, org_category

CLS = {1: "Mainly alpha", 2: "Mainly beta", 3: "Alpha beta", 4: "Few secondary structures"}
DIS = ["mostly disordered", "long disordered region", "neither"]
FOLDS = {"3.40.50": "Rossmann fold", "3.20.20": "TIM barrel", "2.30.30.40": "SH3", "1.10.238.10": "EF-hand",
         "1.10.510.10": "protein kinase C-lobe", "2.30.29.30": "PH domain"}
THRESH = [0.3, 0.5, 0.9, 1.0]


def wilson(k, n, z=1.959964):
    if n == 0:
        return [None, None]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(max(0.0, c - h), 4), round(min(1.0, c + h), 4)]


def prop(k, n):
    return {"k": int(k), "n": int(n), "share": round(k / n, 4) if n else None, "ci95": wilson(k, n)}


def q(x):
    x = np.asarray(x, float)
    return {k: float(np.percentile(x, v)) for k, v in [("p10", 10), ("p25", 25), ("median", 50), ("p75", 75), ("p90", 90)]} \
        | {"n": int(len(x)), "share_gt500": round(float((x > 500).mean()), 4), "max": int(x.max())}


def sfset(s):
    return sorted({x.strip() for x in str(s).split(";") if x.strip() and x.strip() != "nan"})


# ---------------------------------------------------------------- inputs
cols = ["domain", "C", "A", "T", "H", "S35", "S60", "S95", "S100", "S100n", "length", "res"]
cath = pd.read_csv(DL / "cath" / "cath-domain-list.txt", sep=r"\s+", comment="#", names=cols, usecols=range(12))
cath["sf"] = cath.C.astype(str) + "." + cath.A.astype(str) + "." + cath["T"].astype(str) + "." + cath.H.astype(str)
s35 = cath.drop_duplicates(["sf", "S35"])
sf_rank = s35.groupby("sf").size().sort_values(ascending=False)
sf_rank = pd.Series(np.arange(1, len(sf_rank) + 1), index=sf_rank.index)
cs = pd.read_csv(TAB / "cath_s35_metapredict.csv.gz")

mp = pd.read_csv(TAB / "per_protein_mapping.csv", keep_default_na=False, na_values=[""])
mp["clean"] = mp.clean.map({True: True, False: False, "True": True, "False": False})
sm = pd.read_csv(TAB / "saxs_metapredict.csv").set_index("label")
saxs = mp[mp.clean].copy()
saxs = saxs.join(sm, on="label")
ta = pd.read_csv(TAB / "training_annotated.csv.gz", low_memory=False)
counts = json.loads((TAB / "training_counts.json").read_text())
nn = pd.read_csv(TAB / "training_nearest.csv")

sets = {}
sets["SAXS"] = pd.DataFrame({"unit": saxs.uniprot.fillna(saxs.label), "sf": saxs.gene3d_superfamilies.map(sfset),
                             "length": saxs.construct_length, "dis_cat": saxs.dis_cat, "dis_frac": saxs.dis_frac,
                             "org_group": saxs.org_group})
a = ta[ta.component == "afdb"]
sets["AFDB"] = pd.DataFrame({"unit": a.cluster, "sf": a.cath_sf.map(sfset), "length": a.length, "dis_cat": a.dis_cat,
                             "dis_frac": a.dis_frac, "org_group": a.org_group})
m = ta[ta.component.isin(["cath1", "cath2"])]
sets["CATH MD"] = pd.DataFrame({"unit": m.system, "sf": m.cath_sf.map(sfset), "length": m.length, "dis_cat": m.dis_cat,
                                "dis_frac": m.dis_frac, "org_group": m.org_group})
g = ta[ta.component.isin(["megascale", "megasim"])].drop_duplicates("sequence")
sets["MEGAscale"] = pd.DataFrame({"unit": g.system, "sf": g.cath_sf.map(sfset), "length": g.length, "dis_cat": g.dis_cat,
                                  "dis_frac": g.dis_frac, "org_group": g.org_group})
OUT = {"sets": {k: int(len(v)) for k, v in sets.items()}, "counts": counts}
OUT["sets"]["CATH S35"] = int(len(s35))
OUT["components"] = {c: int(n) for c, n in ta.component.value_counts().items()}
OUT["components"]["afdb_member"] = int(pd.read_csv(TAB / "afdb_members.csv.gz", usecols=["member"]).shape[0])
OUT["components"]["octapeptides"] = counts["octapeptides_systems"]
OUT["components"]["megasim_mutants"] = counts["megasim_mutant_systems"]

# ---------------------------------------------------------------- superfamilies and composition
units = {k: v[["unit", "sf"]].explode("sf").dropna().drop_duplicates() for k, v in sets.items()}
units["CATH S35"] = s35[["S35", "sf"]].rename(columns={"S35": "unit"})
for k, u in units.items():
    u["C"] = u.sf.str.split(".").str[0].astype(int)
    u["arch"] = u.sf.str.split(".").str[:2].str.join(".")
    u["topo"] = u.sf.str.split(".").str[:3].str.join(".")
sf_sets = {k: set(u.sf) for k, u in units.items()}
comp = {}
for k, u in units.items():
    u14 = u[u.C.isin(CLS)]
    n = len(u14)
    comp[k] = {"n_units": int(len(u)), "n_units_classes_1_4": int(n),
               "n_with_superfamily": int((sets[k].sf.str.len() > 0).sum()) if k in sets else int(len(s35)),
               "n_superfamilies": len(sf_sets[k]),
               "class": {CLS[c]: prop((u14.C == c).sum(), n) for c in CLS},
               "folds": {f"{t} {nm}": prop((u.topo == t).sum() if t.count(".") == 2 else (u.sf == t).sum(), len(u))
                         for t, nm in FOLDS.items()},
               "rank_beyond_200": prop(int(u.sf.map(sf_rank).gt(200).sum()), len(u)),
               "rank_top_50": prop(int(u.sf.map(sf_rank).le(50).sum()), len(u))}
arch_all = sorted(set(units["CATH S35"].arch.value_counts(normalize=True)[lambda x: x >= 0.01].index),
                  key=lambda a_: [int(v) for v in a_.split(".")])
OUT["architecture"] = {"architectures": arch_all,
                       "share": {k: [round(float((u.arch == a_).mean()), 4) for a_ in arch_all] for k, u in units.items()}}
OUT["composition"] = comp
saxs_sf = sf_sets["SAXS"]
train_main = sf_sets["AFDB"] | sf_sets["CATH MD"] | sf_sets["MEGAscale"]
md = set(units["CATH MD"].sf) | set().union(*ta[ta.component == "megasim"].cath_sf.map(lambda s: set(sfset(s))))
other_sf = set().union(*ta[ta.component == "other_md"].cath_sf.map(lambda s: set(sfset(s))))
train_main |= other_sf
OUT["superfamily_coverage"] = {
    "saxs_superfamilies": len(saxs_sf),
    "saxs_in": {"AFDB": prop(len(saxs_sf & sf_sets["AFDB"]), len(saxs_sf)),
                "CATH MD": prop(len(saxs_sf & sf_sets["CATH MD"]), len(saxs_sf)),
                "MEGAscale": prop(len(saxs_sf & sf_sets["MEGAscale"]), len(saxs_sf)),
                "MD (CATH MD and MegaSim)": prop(len(saxs_sf & md), len(saxs_sf)),
                "any training component": prop(len(saxs_sf & train_main), len(saxs_sf))},
    "saxs_not_in_training": sorted(saxs_sf - train_main),
    "saxs_in_md_only_via_afdb": len((saxs_sf & sf_sets["AFDB"]) - md),
    "training_in_saxs": {k: prop(len(sf_sets[k] & saxs_sf), len(sf_sets[k])) for k in ("AFDB", "CATH MD", "MEGAscale")},
    "cath_superfamilies": int(cath.sf.nunique()),
    "cath_in": {k: prop(len(sf_sets[k] & set(cath.sf)), int(cath.sf.nunique())) for k in ("SAXS", "AFDB", "CATH MD", "MEGAscale")},
    "top200_in": {k: prop(len(sf_sets[k] & set(sf_rank.index[:200])), 200) for k in ("SAXS", "AFDB", "CATH MD", "MEGAscale")},
    "beyond2000_in": {k: prop(len(sf_sets[k] & set(sf_rank.index[2000:])), len(sf_rank) - 2000)
                      for k in ("SAXS", "AFDB", "CATH MD", "MEGAscale")}}
# SAXS entries whose construct carries a superfamily absent from MD training
sx = saxs.assign(sfl=saxs.gene3d_superfamilies.map(sfset))
OUT["superfamily_coverage"]["saxs_entries_with_sf"] = int((sx.sfl.str.len() > 0).sum())
OUT["superfamily_coverage"]["saxs_entries_sf_in_md"] = int(sx.sfl.map(lambda s: bool(set(s) & md)).sum())
OUT["superfamily_coverage"]["saxs_entries_sf_in_any"] = int(sx.sfl.map(lambda s: bool(set(s) & train_main)).sum())

# ---------------------------------------------------------------- length, domains, disorder, organism
OUT["length"] = {k: q(v.length) for k, v in sets.items()} | {"CATH S35": q(cs.length)}
OUT["octapeptides"] = {"systems": counts["octapeptides_systems"], "length": 8}
OUT["multi_superfamily"] = {k: prop(int((v.sf.str.len() >= 2).sum()), len(v)) for k, v in sets.items()}
dis = {}
for k, v in list(sets.items()) + [("CATH S35", cs)]:
    n = len(v)
    dis[k] = {c: prop(int((v.dis_cat == c).sum()), n) for c in DIS} | {"median_dis_frac": round(float(v.dis_frac.median()), 4)}
dis["SAXS_by_class"] = {c: {d_: prop(int((g_.dis_cat == d_).sum()), len(g_)) for d_ in DIS}
                        for c, g_ in saxs.groupby("disorder_class")}
dis["predictor"] = "metapredict v3, residue score >= 0.5; mostly: >= 50% of residues; long region: >= 30 consecutive, under 50%"
OUT["disorder"] = dis
# MEGAscale domains called mostly disordered: natural against de novo designed, and their lengths. The domains were
# selected for folding-stability measurements, so the call is examined by origin and by length
mg = sets["MEGAscale"].assign(designed=sets["MEGAscale"].org_group.eq("Designed"),
                              called=sets["MEGAscale"].dis_cat.eq("mostly disordered"))
tab_ = [[int((mg.designed & mg.called).sum()), int((mg.designed & ~mg.called).sum())],
        [int((~mg.designed & mg.called).sum()), int((~mg.designed & ~mg.called).sum())]]
OUT["megascale_mostly_disordered"] = {
    "natural": prop(int((~mg.designed & mg.called).sum()), int((~mg.designed).sum())),
    "designed": prop(int((mg.designed & mg.called).sum()), int(mg.designed.sum())),
    "fisher_designed_vs_natural_p": float(stats.fisher_exact(tab_)[1]),
    "length_called": q(mg.length[mg.called]) | {"min": int(mg.length[mg.called].min())},
    "length_not_called": q(mg.length[~mg.called]) | {"min": int(mg.length[~mg.called].min())},
    "length_mannwhitney_p": float(stats.mannwhitneyu(mg.length[mg.called], mg.length[~mg.called]).pvalue),
    "length_natural": q(mg.length[~mg.designed]), "length_designed": q(mg.length[mg.designed]),
    "called": [{"system": u, "length": int(n_), "dis_frac": round(float(f_), 3), "designed": bool(d_)}
               for u, n_, f_, d_ in mg.loc[mg.called, ["unit", "length", "dis_frac", "designed"]].itertuples(index=False)]}
pdb = pd.read_csv(TAB / "rcsb_nr30_representatives.csv.gz")
pdb = pdb[pdb.length >= 30]
pdb["org_group"] = [org_category(t, n, (lin if isinstance(lin, str) else "").split("|"))
                    for t, n, lin in zip(pdb.taxid, pdb.sci_name, pdb.lineage_names, strict=True)]
groups = ORG_CATS + ["Designed"]
org = {}
for k, v in list(sets.items()) + [("PDB 30%", pdb)]:
    n = len(v)
    org[k] = {gname: prop(int((v.org_group == gname).sum()), n) for gname in groups}
OUT["organism"] = org

# ---------------------------------------------------------------- overlap (clean entries; all 439 in the table)
c = nn[nn.clean == True]  # noqa: E712
ov = {"n": int(len(c)), "n_all": int(len(nn)), "min_construct_coverage": 0.5}
ov["nearest_identity_quantiles"] = {k: float(np.percentile(c.nearest_identity, v)) for k, v in
                                    [("p10", 10), ("p25", 25), ("median", 50), ("p75", 75), ("p90", 90)]}
ov["at_least"] = {f"{int(t * 100)}": prop(int((c.nearest_identity >= t - 1e-9).sum()), len(c)) for t in THRESH}
ov["at_least_all439"] = {f"{int(t * 100)}": prop(int((nn.nearest_identity >= t - 1e-9).sum()), len(nn)) for t in THRESH}
by = {}
for comp_ in ["afdb", "afdb_member", "cath1", "cath2", "megasim", "megascale", "other_md", "desres"]:
    col = c[f"id_{comp_}"]
    by[comp_] = {f"{int(t * 100)}": prop(int((col >= t - 1e-9).sum()), len(c)) for t in THRESH}
by["md"] = {f"{int(t * 100)}": prop(int((c.nearest_identity_md >= t - 1e-9).sum()), len(c)) for t in THRESH}
ov["by_component"] = by
bins = pd.cut(c.nearest_identity, [-0.01, 0.3, 0.5, 0.9, 1.01], right=False, labels=["<30", "30-50", "50-90", ">=90"])
ov["bins"] = {str(b): int(n) for b, n in bins.value_counts().sort_index().items()}
ov["bins_by_class"] = {cl: {str(b): int(n) for b, n in bins[c.disorder_class == cl].value_counts().sort_index().items()}
                       for cl in ["folded", "partly disordered", "disordered"]}
ov["nearest_component_counts"] = c.nearest_component.value_counts().to_dict()
ov["uniprot_in_afdb_clusters"] = prop(int((c.uniprot_in_afdb_clusters == True).sum()), int(c.uniprot.notna().sum()))  # noqa: E712
ov["uniprot_is_afdb_input"] = prop(int((c.uniprot_is_afdb_input == True).sum()), int(c.uniprot.notna().sum()))  # noqa: E712
ov["sf_group"] = c.sf_group.value_counts().to_dict()
ov["sf_group_md"] = c.sf_group_md.value_counts().to_dict()
near = nn[nn.nearest_identity >= 0.9].sort_values("nearest_identity", ascending=False)
ov["near_identical"] = [{"label": r["label"], "clean": bool(r["clean"]), "identity": round(float(r["nearest_identity"]), 3),
                         "component": r["nearest_component"], "hit": r["nearest_hit"],
                         "coverage": round(float(r[f"cov_{r['nearest_component']}"]), 3),
                         "construct_length": int(r["construct_length"]), "disorder_class": r["disorder_class"],
                         "uniprot": r["uniprot"] if isinstance(r["uniprot"], str) else None}
                        for _, r in near.iterrows()]
OUT["overlap"] = ov
(TAB / "training_summary.json").write_text(json.dumps(OUT, indent=1, default=str))
print(json.dumps({k: OUT[k] for k in ("sets", "superfamily_coverage")}, indent=1, default=str)[:3000])
print(json.dumps(ov["at_least"], indent=1), ov["bins"], ov["bins_by_class"])
