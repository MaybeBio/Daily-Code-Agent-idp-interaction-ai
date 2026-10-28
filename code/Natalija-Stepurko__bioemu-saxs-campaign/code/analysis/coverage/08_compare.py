"""All eight comparisons. Writes tables/summary.json plus comparison tables used by 09_figures.py.

Units
- Ours: the 399 clean entries (bioemu, clean == True). Domain-level comparisons use unique
  (UniProt accession, CATH superfamily) pairs from CATH-Gene3D hits inside the construct
  region, the protein-family analogue of a CATH S35 sequence cluster.
- CATH: v4.4.0 domain list; class/architecture/superfamily composition counted over S35
  sequence clusters (one representative per cluster).
- Intervals: Wilson 95% for single proportions; cluster bootstrap over UniProt accessions
  (2000 resamples) for composition shares and enrichment ratios.
"""
import json, re
import numpy as np
import pandas as pd
from scipy import stats
from common import *

rng = np.random.default_rng(0)
B = 2000
OUT = {}


def wilson(k, n, z=1.96):
    if n == 0:
        return [np.nan, np.nan]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 4), round(c + h, 4)]


def prop(k, n):
    return {"k": int(k), "n": int(n), "share": round(k / n, 4) if n else None, "ci95": wilson(k, n)}


# ---------------------------------------------------------------- CATH reference
def read_domain_list(path):
    cols = ["domain", "C", "A", "T", "H", "S35", "S60", "S95", "S100", "S100n", "length", "res"]
    d = pd.read_csv(path, sep=r"\s+", comment="#", names=cols, usecols=range(12))
    d["arch"] = d.C.astype(str) + "." + d.A.astype(str)
    d["topo"] = d.arch + "." + d["T"].astype(str)
    d["sf"] = d.topo + "." + d["H"].astype(str)
    return d


cath = read_domain_list(DL / "cath" / "cath-domain-list.txt")
cath_ver = re.search(r"CATH VERSION:\s*(\S+)", (DL / "cath" / "cath-domain-list.txt").read_text()[:2000]).group(1)
s35 = cath.drop_duplicates(["sf", "S35"])  # one row per S35 cluster
cath43 = read_domain_list(DL / "cath" / "cath-domain-list-S35-v4_3_0.txt")
record_version("CATH", {"version": cath_ver, "version_date": "16.12.2024",
                        "source": "https://download.cathdb.info/cath/releases/latest-release/cath-classification-data/",
                        "n_domains": len(cath), "n_s35_clusters": len(s35)})
names = {}
for line in (DL / "cath" / "cath-names.txt").read_text().splitlines():
    if line.startswith("#") or not line.strip():
        continue
    node, rest = line.split(None, 1)
    names[node] = rest.split(":", 1)[1].strip() if ":" in rest else ""

sf_s35 = s35.groupby("sf").size().sort_values(ascending=False)
sf_rank = pd.Series(np.arange(1, len(sf_s35) + 1), index=sf_s35.index)
OUT["cath_reference"] = {"version": cath_ver, "n_domains": len(cath), "n_s35_clusters": len(s35),
                         "n_superfamilies": cath.sf.nunique(), "n_topologies": cath.topo.nunique(),
                         "n_architectures": cath.arch.nunique(), "n_classes": int(cath.C.nunique()),
                         "v4_3_0": {"n_s35_clusters": len(cath43), "n_superfamilies": cath43.sf.nunique(),
                                    "n_topologies": cath43.topo.nunique(), "n_architectures": cath43.arch.nunique()}}

# ---------------------------------------------------------------- the SAXS-set table
df = pd.read_csv(TAB / "per_protein_mapping.csv", keep_default_na=False, na_values=[""])
for c in ["gene3d_superfamilies", "ted_superfamilies", "pfam"]:
    df[c] = df[c].fillna("")
for c in ["has_pdb_50", "has_pdb_any_overlap", "clean", "has_cath_gene3d", "has_cath_ted", "has_cath_any"]:
    df[c] = df[c].map({True: True, False: False, "True": True, "False": False})
C = df[df.clean].copy()
N = len(C)
OUT["n_entries_all"] = len(df); OUT["n_clean"] = N

# ---------------------------------------------------------------- mapping success
mp = {}
for name, sub in [("clean", C), ("all", df)]:
    mp[name] = {"n": len(sub), "mapped": prop((sub.map_method != "unmapped").sum(), len(sub)),
                "by_method": sub.map_method.value_counts().to_dict(),
                "sasbdb_accession_replaced": int(sub.map_detail.fillna("").str.contains("SASBDB accession").sum()),
                "distinct_accessions": int(sub.uniprot.nunique()),
                "median_alignment_identity": round(float(sub.aln_identity.median()), 4),
                "n_identity_below_0.9": int((sub.aln_identity < 0.9).sum()),
                "n_construct_coverage_below_0.8": int((sub.construct_coverage < 0.8).sum())}
mp["clean"]["unmapped_labels"] = C[C.map_method == "unmapped"][["label", "molecule_name"]].values.tolist()
OUT["mapping"] = mp

# ---------------------------------------------------------------- 1 coverage
def sfs(s):
    return [x for x in s.split(";") if x]


C["g3_list"] = C.gene3d_superfamilies.map(sfs)
C["ted_list"] = C.ted_superfamilies.map(sfs)
ours_sf = sorted({x for l in C.g3_list for x in l})
ours_sf_ted = sorted({x for l in C.ted_list for x in l})
ours_sf_union = sorted(set(ours_sf) | set(ours_sf_ted))
cath_sfs = set(cath.sf)


def levels(sfl):
    return ({".".join(x.split(".")[:3]) for x in sfl}, {".".join(x.split(".")[:2]) for x in sfl})


t_g3, a_g3 = levels(ours_sf); t_u, a_u = levels(ours_sf_union)
top = {k: set(sf_s35.index[:k]) for k in (50, 200)}
cov = {"has_gene3d_superfamily": prop(C.has_cath_gene3d.sum(), N),
       "has_ted_superfamily": prop(C.has_cath_ted.sum(), N),
       "has_gene3d_or_ted": prop(C.has_cath_any.sum(), N),
       "has_gene3d_by_disorder_class": {k: prop(g.has_cath_gene3d.sum(), len(g)) for k, g in C.groupby("disorder_class")},
       "gene3d": {"distinct_superfamilies": len(ours_sf), "distinct_topologies": len(t_g3),
                  "distinct_architectures": len(a_g3),
                  "superfamilies_not_in_cath_v4.4": sorted(set(ours_sf) - cath_sfs),
                  "share_of_cath_superfamilies": round(len(set(ours_sf) & cath_sfs) / cath.sf.nunique(), 4),
                  "share_of_cath_topologies": round(len(t_g3 & set(cath.topo)) / cath.topo.nunique(), 4),
                  "share_of_cath_architectures": round(len(a_g3 & set(cath.arch)) / cath.arch.nunique(), 4),
                  "share_of_cath_s35_clusters_in_represented_superfamilies":
                      round(float(sf_s35[sf_s35.index.isin(ours_sf)].sum() / sf_s35.sum()), 4)},
       "gene3d_or_ted": {"distinct_superfamilies": len(ours_sf_union), "distinct_topologies": len(t_u),
                         "distinct_architectures": len(a_u),
                         "share_of_cath_s35_clusters_in_represented_superfamilies":
                             round(float(sf_s35[sf_s35.index.isin(ours_sf_union)].sum() / sf_s35.sum()), 4)}}
for k in (50, 200):
    cov[f"top{k}_superfamilies_present_gene3d"] = prop(len(top[k] & set(ours_sf)), k)
    cov[f"top{k}_superfamilies_present_gene3d_or_ted"] = prop(len(top[k] & set(ours_sf_union)), k)
OUT["coverage"] = cov
pd.DataFrame({"sf": sf_s35.index, "s35_clusters": sf_s35.values, "rank": sf_rank.values,
              "in_ours_gene3d": sf_s35.index.isin(ours_sf), "in_ours_ted": sf_s35.index.isin(ours_sf_ted),
              "name": [names.get(x, "") for x in sf_s35.index]}).to_csv(TAB / "cath_superfamily_rank.csv", index=False)

# ---------------------------------------------------------------- 2 composition
units = (C[["uniprot", "g3_list"]].explode("g3_list").dropna().drop_duplicates()
         .rename(columns={"g3_list": "sf"}))
units["C"] = units.sf.str.split(".").str[0].astype(int)
units["arch"] = units.sf.str.split(".").str[:2].str.join(".")
units["topo"] = units.sf.str.split(".").str[:3].str.join(".")
n_inst = int(C.n_gene3d_domains.sum())
CLS = {1: "Mainly alpha", 2: "Mainly beta", 3: "Alpha beta", 4: "Few secondary structures"}
s35c = s35[s35.C.isin(CLS)]
ref_cls = s35c.C.value_counts(normalize=True)
u14 = units[units.C.isin(CLS)]
accs = u14.uniprot.unique()
by_acc = {a: g for a, g in u14.groupby("uniprot")}


def boot_shares(col, cats):
    obs = u14[col].value_counts(normalize=True).reindex(cats, fill_value=0)
    bs = []
    for _ in range(B):
        samp = pd.concat([by_acc[a] for a in rng.choice(accs, len(accs))])
        bs.append(samp[col].value_counts(normalize=True).reindex(cats, fill_value=0).values)
    bs = np.array(bs)
    return obs, np.percentile(bs, 2.5, 0), np.percentile(bs, 97.5, 0)


obs, lo, hi = boot_shares("C", list(CLS))
cls_tab = pd.DataFrame({"class": list(CLS), "name": list(CLS.values()), "ours_units": u14.C.value_counts().reindex(list(CLS), fill_value=0).values,
                        "ours_share": obs.values, "ours_lo": lo, "ours_hi": hi,
                        "cath_s35": s35c.C.value_counts().reindex(list(CLS)).values,
                        "cath_share": ref_cls.reindex(list(CLS)).values})
cls_tab["enrichment"] = cls_tab.ours_share / cls_tab.cath_share
cls_tab["enr_lo"] = cls_tab.ours_lo / cls_tab.cath_share
cls_tab["enr_hi"] = cls_tab.ours_hi / cls_tab.cath_share
cls_tab.to_csv(TAB / "class_composition.csv", index=False)

arch_ref = s35c.arch.value_counts(normalize=True)
arch_keep = sorted(set(arch_ref[arch_ref >= 0.01].index) | set(u14.arch.value_counts()[lambda x: x >= 3].index),
                   key=lambda a: [int(v) for v in a.split(".")])
obs, lo, hi = boot_shares("arch", arch_keep)
arch_tab = pd.DataFrame({"arch": arch_keep, "name": [names.get(a, "") for a in arch_keep],
                         "ours_units": u14.arch.value_counts().reindex(arch_keep, fill_value=0).values,
                         "ours_share": obs.values, "ours_lo": lo, "ours_hi": hi,
                         "cath_s35": s35c.arch.value_counts().reindex(arch_keep, fill_value=0).values,
                         "cath_share": arch_ref.reindex(arch_keep, fill_value=0).values})
arch_tab["enrichment"] = arch_tab.ours_share / arch_tab.cath_share
arch_tab["enr_lo"] = arch_tab.ours_lo / arch_tab.cath_share
arch_tab["enr_hi"] = arch_tab.ours_hi / arch_tab.cath_share
arch_tab.to_csv(TAB / "architecture_composition.csv", index=False)
OUT["composition"] = {"unit": "unique (UniProt accession, Gene3D superfamily) pairs inside constructs, classes 1-4",
                      "n_units": int(len(u14)), "n_units_all_classes": int(len(units)),
                      "n_domain_instances": n_inst, "n_proteins": int(len(accs)),
                      "cath_s35_clusters_classes_1_4": int(len(s35c)),
                      "class": cls_tab.round(4).to_dict("records"),
                      "architecture": arch_tab.round(4).to_dict("records")}

# ---------------------------------------------------------------- 3 over/under superfamilies and topologies
def over_under(level, ref_counts):
    obs = units.groupby(level).uniprot.nunique()
    n = len(units)
    exp = ref_counts / ref_counts.sum() * n
    t = pd.DataFrame({"observed": obs.reindex(exp.index, fill_value=0), "expected": exp,
                      "cath_s35": ref_counts})
    extra = obs.index.difference(exp.index)
    t["ratio"] = t.observed / t.expected
    t["p_over"] = stats.poisson.sf(t.observed - 1, t.expected)
    t["p_under"] = stats.poisson.cdf(t.observed, t.expected)
    t["name"] = [names.get(x, "") for x in t.index]
    over = t[t.observed >= 3].sort_values("p_over").head(15)
    under = t.sort_values("p_under").head(15)
    return t, over, under, list(extra)


sf_t, sf_over, sf_under, sf_extra = over_under("sf", sf_s35)
tp_t, tp_over, tp_under, _ = over_under("topo", s35.groupby("topo").size())
sf_t.to_csv(TAB / "superfamily_observed_expected.csv")
tp_t.to_csv(TAB / "topology_observed_expected.csv")
fmt = lambda t: [{"id": i, "name": r["name"], "observed": int(r.observed), "expected": round(r.expected, 2),
                  "cath_s35": int(r.cath_s35), "p": float(f"{min(r.p_over, r.p_under):.2g}")} for i, r in t.iterrows()]
OUT["over_under"] = {"n_units": int(len(units)), "superfamily_over": fmt(sf_over.head(10)),
                     "superfamily_under": fmt(sf_under.head(10)), "topology_over": fmt(tp_over.head(10)),
                     "topology_under": fmt(tp_under.head(10)), "ours_sf_absent_from_cath_v4.4": sf_extra,
                     "note": "expected = n_units x CATH S35 share; Poisson tail p, uncorrected, for ranking only"}

# ---------------------------------------------------------------- 4 organism
sp = pd.read_csv(TAB / "swissprot_summary.csv.gz")
pdb = pd.read_csv(TAB / "rcsb_nr30_representatives.csv.gz")
pdb["org_group"] = [org_category(t, n, (l if isinstance(l, str) else "").split("|"))
                    for t, n, l in zip(pdb.taxid, pdb.sci_name, pdb.lineage_names)]
pdb30 = pdb[pdb.length >= 30]  # drop short peptides (20.6% of experimental NR30 representatives)
org = pd.DataFrame({"ours": C.org_group.value_counts(), "pdb_nr30": pdb30.org_group.value_counts(),
                    "swissprot": sp.org_group.value_counts()}).reindex(ORG_CATS).fillna(0).astype(int)
orgs = org / org.sum()
org_tab = orgs.add_suffix("_share").join(org)
org_tab["ours_lo"] = [wilson(k, N)[0] for k in org.ours]; org_tab["ours_hi"] = [wilson(k, N)[1] for k in org.ours]
org_tab.to_csv(TAB / "organism_comparison.csv")
OUT["organism"] = {"groups": org_tab.round(4).reset_index(names="group").to_dict("records"),
                   "pdb_nr30_ge30_no_source_organism": int(pdb30.taxid.isna().sum()),
                   "n_pdb_nr30_experimental": len(pdb), "n_pdb_nr30_ge30": len(pdb30), "n_swissprot": len(sp),
                   "ours_other_top": C[C.org_group == "Other"].organism.value_counts().head(10).to_dict()}

# ---------------------------------------------------------------- 5 length, domains per protein
q = lambda x: {k: float(np.percentile(x, v)) for k, v in [("p10", 10), ("p25", 25), ("median", 50), ("p75", 75), ("p90", 90)]}
pdbL = pdb.length.dropna()
OUT["length"] = {"ours_construct": q(C.construct_length), "ours_fulllength_parent": q(C.uniprot_length.dropna()),
                 "pdb_nr30_all": q(pdbL), "pdb_nr30_ge30": q(pdbL[pdbL >= 30]), "swissprot": q(sp.length),
                 "pdb_nr30_share_lt30": round(float((pdbL < 30).mean()), 4),
                 "ks_ours_vs_pdb_nr30_ge30": float(f"{stats.ks_2samp(C.construct_length, pdbL[pdbL >= 30]).pvalue:.3g}"),
                 "ks_ours_vs_swissprot": float(f"{stats.ks_2samp(C.construct_length, sp.length).pvalue:.3g}"),
                 "ours_share_gt500": prop((C.construct_length > 500).sum(), N),
                 "pdb_nr30_ge30_share_gt500": round(float((pdbL[pdbL >= 30] > 500).mean()), 4),
                 "swissprot_share_gt500": round(float((sp.length > 500).mean()), 4)}
Cm = C[C.map_method != "unmapped"]
bins = lambda x: pd.cut(x, [-1, 0, 1, 2, 3, 1e9], labels=["0", "1", "2", "3", "4+"]).value_counts(normalize=True).sort_index().round(4).to_dict()
OUT["domains_per_protein"] = {"ours_construct_gene3d": bins(Cm.n_gene3d_domains), "ours_construct_ted": bins(Cm.n_ted_domains),
                              "ours_parent_fulllength_gene3d": bins(Cm.n_gene3d_domains_fulllength),
                              "swissprot_fulllength_gene3d": bins(sp.n_gene3d_matches),
                              "ours_construct_gene3d_mean": round(float(Cm.n_gene3d_domains.mean()), 3),
                              "ours_parent_fulllength_gene3d_mean": round(float(Cm.n_gene3d_domains_fulllength.mean()), 3),
                              "swissprot_fulllength_gene3d_mean": round(float(sp.n_gene3d_matches.mean()), 3)}
dpp = pd.DataFrame({"ours_construct": bins(Cm.n_gene3d_domains), "ours_parent_fulllength": bins(Cm.n_gene3d_domains_fulllength),
                    "swissprot": bins(sp.n_gene3d_matches)})
dpp.to_csv(TAB / "domains_per_protein.csv")

# ---------------------------------------------------------------- 6 disorder
def dis_cat(frac, longest):
    return np.where(frac >= 0.5, "mostly disordered (>=50% residues)",
                    np.where(longest >= 30, "long IDR (>=30 aa), <50%", "no long IDR"))


sp["dcat"] = dis_cat(sp.dis_res / sp.length, sp.longest_dis)
Cm = Cm.assign(dcat_construct=dis_cat(Cm.up_disordered_frac_construct, Cm.up_longest_disordered_construct),
               dcat_full=dis_cat(Cm.up_disordered_frac_fulllength, Cm.up_longest_disordered_fulllength))
DC = ["mostly disordered (>=50% residues)", "long IDR (>=30 aa), <50%", "no long IDR"]
dis = pd.DataFrame({"ours_construct": Cm.dcat_construct.value_counts(), "ours_parent_fulllength": Cm.dcat_full.value_counts(),
                    "swissprot_all": sp.dcat.value_counts(), "swissprot_with_pdb": sp[sp.n_pdb > 0].dcat.value_counts(),
                    "swissprot_without_pdb": sp[sp.n_pdb == 0].dcat.value_counts()}).reindex(DC).fillna(0).astype(int)
dis.to_csv(TAB / "disorder_comparison.csv")
nm = len(Cm)
OUT["disorder"] = {
    "ours_adopt2_class_clean": {k: prop(v, N) for k, v in C.disorder_class.value_counts().items()},
    "uniprot_disordered_annotation": {
        "definition": "UniProt REGION features with note 'Disordered' (MobiDB-lite predictions plus curated); "
                      "mostly = >=50% of residues; long IDR = >=30 consecutive residues",
        "ours_construct": {k: prop(dis.loc[k, "ours_construct"], nm) for k in DC},
        "ours_parent_fulllength": {k: prop(dis.loc[k, "ours_parent_fulllength"], nm) for k in DC},
        "swissprot_all": {k: round(dis.loc[k, "swissprot_all"] / len(sp), 4) for k in DC},
        "swissprot_with_pdb_xref": {k: round(dis.loc[k, "swissprot_with_pdb"] / (sp.n_pdb > 0).sum(), 4) for k in DC},
        "swissprot_without_pdb_xref": {k: round(dis.loc[k, "swissprot_without_pdb"] / (sp.n_pdb == 0).sum(), 4) for k in DC},
        "n_swissprot_with_pdb": int((sp.n_pdb > 0).sum())},
    "crosstab_adopt2_vs_uniprot_construct": pd.crosstab(Cm.disorder_class, Cm.dcat_construct).to_dict()}

# ---------------------------------------------------------------- 7 structure availability
st = C[C.map_method != "unmapped"].copy()
for c_ in ["has_pdb_50", "has_pdb_any_overlap", "uniprot_has_pdb_xref"]:
    st[c_] = st[c_].map({True: True, False: False, "True": True, "False": False}).astype(bool)
ct = pd.crosstab(st.disorder_class, st.has_pdb_50)
dd = st.disorder_class != "folded"
fis = stats.fisher_exact([[int((~st.has_pdb_50 & dd).sum()), int((st.has_pdb_50 & dd).sum())],
                          [int((~st.has_pdb_50 & ~dd).sum()), int((st.has_pdb_50 & ~dd).sum())]])
nostruct = st[~st.has_pdb_50]
OUT["structure_availability"] = {
    "definition": "SIFTS best_structures chains from atomic-model methods (X-ray, EM, NMR, ...; SAXS models excluded); "
                  "UniProt residue ranges are SEQRES-level and can include unmodelled residues",
    "n_mapped_clean": len(st),
    "no_pdb_ge50pct_of_construct": prop((~st.has_pdb_50).sum(), len(st)),
    "no_pdb_overlap_at_all": prop((~st.has_pdb_any_overlap).sum(), len(st)),
    "uniprot_entry_without_any_pdb_xref": prop((~st.uniprot_has_pdb_xref.astype(bool)).sum(), len(st)),
    "by_disorder_class_no_pdb_ge50": {k: prop((~g.has_pdb_50).sum(), len(g)) for k, g in st.groupby("disorder_class")},
    "composition_of_no_structure_set": nostruct.disorder_class.value_counts().to_dict(),
    "fisher_disordered_or_partly_vs_folded": {"odds_ratio": round(fis[0], 3), "p": float(f"{fis[1]:.3g}")},
    "unmapped_clean_entries_not_assessed": int((C.map_method == "unmapped").sum()),
    "swissprot_share_with_pdb_xref": round(float((sp.n_pdb > 0).mean()), 4)}
ct.to_csv(TAB / "structure_by_disorder.csv")

# ---------------------------------------------------------------- 8 redundancy
acc_counts = st.uniprot.value_counts()
sf_entries = C.explode("g3_list").dropna(subset=["g3_list"]).groupby("g3_list").label.nunique().sort_values(ascending=False)
sf_accs = C.explode("g3_list").dropna(subset=["g3_list"]).groupby("g3_list").uniprot.nunique()
e_sf = C[C.has_cath_gene3d].copy()
share_other = []
for _, r in e_sf.iterrows():
    others = e_sf[(e_sf.uniprot != r.uniprot)]
    share_other.append(any(set(r.g3_list) & set(l) for l in others.g3_list))
OUT["redundancy"] = {
    "distinct_accessions_mapped_clean": int(st.uniprot.nunique()), "mapped_clean_entries": len(st),
    "entries_sharing_accession_with_another_entry": prop(int(acc_counts[acc_counts > 1].sum()), len(st)),
    "accessions_with_multiple_entries": int((acc_counts > 1).sum()),
    "top_accessions": [{"uniprot": a, "entries": int(n), "name": str(st[st.uniprot == a].molecule_name.iloc[0])}
                       for a, n in acc_counts.head(10).items()],
    "entries_with_gene3d": int(len(e_sf)),
    "entries_sharing_a_superfamily_with_an_entry_of_another_accession": prop(int(sum(share_other)), len(e_sf)),
    "top_superfamilies_by_entries": [{"sf": s, "name": names.get(s, ""), "entries": int(n), "accessions": int(sf_accs[s])}
                                     for s, n in sf_entries.head(10).items()],
    "superfamilies_seen_once": int((sf_entries == 1).sum())}

# Pfam
pf = C.pfam.map(sfs)
OUT["pfam"] = {"entries_with_pfam_in_construct": prop(int((pf.str.len() > 0).sum()), N),
               "distinct_pfam_families": len({x for l in pf for x in l})}

(TAB / "summary.json").write_text(json.dumps({"versions": json.loads(VERSIONS.read_text()), **OUT}, indent=1, default=str))
print(json.dumps(OUT, indent=1, default=str)[:6000])
