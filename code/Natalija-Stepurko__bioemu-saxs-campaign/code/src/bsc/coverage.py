"""Stage `coverage`: where the SAXS set sits in known protein structure space, and whether BioEmu-1 fits
profiles better when the construct has an experimental structure.

Reads only tracked outputs: results/coverage/ (written by analysis/coverage, which needs network access
and its own requirements) and results/resolvability.csv. Writes results/coverage/structure_split.json,
results/coverage/rank_coverage.json and the coverage figures, and, from the training-data outputs
(training_nearest.csv, training_summary.json), results/coverage/training_closeness.json (BioEmu-1's raw
error against the closeness of each construct to the training data), results/coverage/category_evidence.json
(the error on candidate protein categories for data collection, and the share of all archive entries that pass
the quality checks by disorder class, from results/entries.csv) and three training-data figures.

Structure groups, from SIFTS residue ranges of atomic-model PDB chains mapped onto the construct region
of the UniProt entry (analysis/coverage/07_build_protein_table.py):
  pdb50    experimental chains cover at least 50% of the construct (has_pdb_50)
  partial  some overlap, under 50%
  none     no overlap
Entries without a UniProt mapping are left out. Tests compare pdb50 against the other two (under 50%).
"""
from __future__ import annotations

import json

import matplotlib
import matplotlib.patches
import numpy as np
import pandas as pd
from scipy import stats

from bsc import analyse, report, robust
from bsc import config as C

COV = C.RESULTS / "coverage"
GROUPS = ["pdb50", "partial", "none"]
GROUP_LABEL = {"pdb50": "PDB ≥ 50%", "partial": "partial", "none": "no PDB"}
RAW_FIT = 2.0
SEED = 20261014
RANK_EDGES = [0, 10, 50, 200, 500, 1000, 2000]


def _bool(s: pd.Series) -> pd.Series:
    return s.map({True: True, False: False, "True": True, "False": False})


def structure_group(mp: pd.DataFrame) -> pd.Series:
    """pdb50 / partial / none per label; NaN for entries without a UniProt mapping."""
    p50, anyo = _bool(mp["has_pdb_50"]), _bool(mp["has_pdb_any_overlap"])
    g = pd.Series(np.where(p50 == True, "pdb50", np.where(anyo == True, "partial", "none")),  # noqa: E712
                  index=mp.index, dtype=object)
    g[(mp["map_method"] == "unmapped") | p50.isna()] = np.nan
    return g


def wilson(k: int, n: int, z: float = 1.959964) -> list[float]:
    """Wilson score interval of a proportion."""
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [float(max(0.0, c - h)), float(min(1.0, c + h))]


def van_elteren(x: np.ndarray, group: np.ndarray, strata: np.ndarray) -> dict:
    """Stratified Wilcoxon rank-sum test (van Elteren): rank sums of group True within each stratum,
    weighted by 1/(n_s + 1), normal approximation with the tie-corrected variance. Two-sided p."""
    t = e = v = 0.0
    for s in np.unique(strata):
        k = strata == s
        xs, gs = x[k], group[k]
        m, n = int(gs.sum()), int(len(gs))
        if m == 0 or m == n:
            continue
        r = stats.rankdata(xs)
        w = 1.0 / (n + 1)
        t += w * r[gs].sum()
        e += w * m * (n + 1) / 2
        _, cnt = np.unique(xs, return_counts=True)
        tie = 1 - (cnt ** 3 - cnt).sum() / (n ** 3 - n)
        v += w * w * m * (n - m) * (n + 1) / 12 * tie
    z = (t - e) / np.sqrt(v) if v > 0 else float("nan")
    return {"z": float(z), "p_two_sided": float(2 * stats.norm.sf(abs(z)))}


def holm(p: dict[str, float]) -> dict[str, float]:
    """Holm step-down adjusted p values."""
    keys = sorted(p, key=p.get)
    out, run = {}, 0.0
    for i, k in enumerate(keys):
        run = max(run, min(1.0, (len(keys) - i) * p[k]))
        out[k] = run
    return out


def summarise(g: pd.DataFrame, seed: int) -> dict:
    """n, raw-fit share (Wilson interval), medians of raw chi2, NRMSD and Rg ratio (bootstrap intervals),
    share whose target fit is not reached (Wilson interval)."""
    n = int(len(g))
    k_fit = int((g["chi2_raw"] <= RAW_FIT).sum())
    k_nf = int((g["resolvability"] == analyse.KIND_NOTFIT).sum())

    def med(col, off):
        v = g[col].dropna().to_numpy(float)
        if len(v) == 0:
            return None, [None, None]
        return float(np.median(v)), list(robust.bootstrap_ci(v, np.median, analyse.N_BOOT, seed + off))

    out = {"n": n, "raw_fit": {"k": k_fit, "share": k_fit / n if n else None, "ci95": wilson(k_fit, n)},
           "target_not_reached": {"k": k_nf, "share": k_nf / n if n else None, "ci95": wilson(k_nf, n)}}
    for col, key, off in (("chi2_raw", "chi2_raw_median", 1), ("nrmsd_raw", "nrmsd_raw_median", 2),
                          ("rg_ratio", "rg_ratio_median", 3)):
        m, c = med(col, off)
        out[key], out[key.replace("_median", "_ci95")] = m, c
    return out


def _mw(a: pd.Series, b: pd.Series) -> float:
    a, b = a.dropna(), b.dropna()
    if not (len(a) and len(b)):
        return float("nan")
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def structure_split(d: pd.DataFrame, mp: pd.DataFrame) -> dict:
    """Raw-fit statistics of BioEmu-1 by structure group, pooled and within each disorder class, with
    two-sided Mann-Whitney tests (pdb50 against under 50%) on raw chi2 and on the Rg ratio, Holm-adjusted
    within the six within-class tests, and the class-stratified van Elteren test."""
    m = d[(d["model"] == C.MODEL) & d["clean"]].copy()
    grp = structure_group(mp.set_index("label"))
    m["structure"] = m["label"].map(grp)
    n_unmapped = int(m["structure"].isna().sum())
    m = m.dropna(subset=["structure"])
    out = {"definition": "SIFTS residue ranges of PDB chains from atomic-model methods mapped onto the "
                         "UniProt region of the construct; pdb50: at least 50% of the construct covered; "
                         "partial: some overlap, under 50%; none: no overlap. Raw fit: reduced chi2 <= 2. "
                         "Share intervals: Wilson 95%; median intervals: percentile bootstrap over entries "
                         "(2,000 resamples).",
           "n": int(len(m)), "n_unmapped_excluded": n_unmapped, "groups": GROUPS}
    seed = SEED
    out["pooled"] = {}
    for gname in GROUPS + ["under50"]:
        g = m[m["structure"] != "pdb50"] if gname == "under50" else m[m["structure"] == gname]
        out["pooled"][gname] = summarise(g, seed)
        seed += 10
    w, wo = m[m["structure"] == "pdb50"], m[m["structure"] != "pdb50"]
    out["pooled_tests"] = {
        "chi2_raw_mannwhitney_p": _mw(w["chi2_raw"], wo["chi2_raw"]),
        "rg_ratio_mannwhitney_p": _mw(w["rg_ratio"], wo["rg_ratio"]),
        "chi2_raw_pdb50_vs_none_p": _mw(w["chi2_raw"], m.loc[m["structure"] == "none", "chi2_raw"]),
        "nrmsd_raw_mannwhitney_p": _mw(w["nrmsd_raw"], wo["nrmsd_raw"])}
    k = m["rg_ratio"].notna().to_numpy()
    out["stratified_by_class"] = {
        "chi2_raw": van_elteren(np.log(m["chi2_raw"].to_numpy(float)), (m["structure"] == "pdb50").to_numpy(),
                                m["disorder_class"].to_numpy()),
        "rg_ratio": van_elteren(m["rg_ratio"].to_numpy(float)[k], (m["structure"] == "pdb50").to_numpy()[k],
                                m["disorder_class"].to_numpy()[k])}
    # composition of the groups by disorder class (the confounder)
    out["class_share_by_group"] = {
        gname: {c: float((m.loc[m["structure"] == gname, "disorder_class"] == c).mean())
                for c in analyse.CLASS_ORDER}
        for gname in GROUPS}
    out["pdb50_share_by_class"] = {
        c: {"k": int((g["structure"] == "pdb50").sum()), "n": int(len(g)),
            "share": float((g["structure"] == "pdb50").mean()),
            "ci95": wilson(int((g["structure"] == "pdb50").sum()), int(len(g)))}
        for c, g in m.groupby("disorder_class")}
    out["by_class"], raw_p = {}, {}
    for cls in analyse.CLASS_ORDER:
        g = m[m["disorder_class"] == cls]
        rec = {}
        for gname in GROUPS + ["under50"]:
            sub = g[g["structure"] != "pdb50"] if gname == "under50" else g[g["structure"] == gname]
            rec[gname] = summarise(sub, seed)
            seed += 10
        a, b = g[g["structure"] == "pdb50"], g[g["structure"] != "pdb50"]
        rec["tests"] = {"chi2_raw_mannwhitney_p": _mw(a["chi2_raw"], b["chi2_raw"]),
                        "rg_ratio_mannwhitney_p": _mw(a["rg_ratio"], b["rg_ratio"]),
                        "nrmsd_raw_mannwhitney_p": _mw(a["nrmsd_raw"], b["nrmsd_raw"])}
        raw_p[f"{cls}|chi2_raw"] = rec["tests"]["chi2_raw_mannwhitney_p"]
        raw_p[f"{cls}|rg_ratio"] = rec["tests"]["rg_ratio_mannwhitney_p"]
        out["by_class"][cls] = rec
    adj = holm(raw_p)
    for key, p in adj.items():
        cls, metric = key.split("|")
        out["by_class"][cls]["tests"][f"{metric}_holm_p"] = p
    out["holm_family"] = "the six within-class tests (three classes, raw chi2 and Rg ratio)"
    return out


def rank_coverage(rk: pd.DataFrame) -> dict:
    """Share of CATH superfamilies present in the SAXS set (CATH-Gene3D within the construct) by rank bin
    of S35 cluster count."""
    edges = RANK_EDGES + [int(rk["rank"].max())]
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        sub = rk[(rk["rank"] > lo) & (rk["rank"] <= hi)]
        k = int(sub["in_ours_gene3d"].sum())
        bins.append({"first": lo + 1, "last": hi, "n": int(len(sub)), "k": k, "share": k / len(sub),
                     "ci95": wilson(k, len(sub)), "s35_clusters": int(sub["s35_clusters"].sum())})
    return {"bins": bins, "n_superfamilies": int(len(rk)),
            "n_present": int(rk["in_ours_gene3d"].sum())}


ID_BINS = [0.0, 0.3, 0.5, 0.9, 1.0001]
ID_LABELS = ["<30", "30-50", "50-90", ">=90"]
CLOSE_METRICS = ("chi2_raw", "raw_fit", "nrmsd_raw", "rg_ratio")
MW_METRICS = ("chi2_raw", "nrmsd_raw", "rg_ratio")


def identity_bin(x: pd.Series) -> pd.Series:
    """Nearest-training identity (0-1) into the four bins <30%, 30-50%, 50-90%, >=90%."""
    return pd.cut(x, ID_BINS, right=False, labels=ID_LABELS).astype(object)


def _spearman(x: pd.Series, y: pd.Series) -> dict:
    k = x.notna() & y.notna()
    if k.sum() < 5 or x[k].nunique() < 2 or y[k].nunique() < 2:
        return {"rho": None, "p": None, "n": int(k.sum())}
    r = stats.spearmanr(x[k], y[k])
    return {"rho": float(r.statistic), "p": float(r.pvalue), "n": int(k.sum())}


def _bioemu_with_training(d: pd.DataFrame, nn: pd.DataFrame) -> pd.DataFrame:
    m = d[(d["model"] == C.MODEL) & d["clean"]].copy()
    t = nn.set_index("label")
    for col in ("nearest_identity", "nearest_identity_md", "sf_group", "sf_group_md", "id_afdb_member"):
        m[col] = m["label"].map(t[col])
    m["identity_bin"] = identity_bin(m["nearest_identity"])
    m["raw_fit"] = (m["chi2_raw"] <= RAW_FIT).astype(float)
    return m


def training_closeness(d: pd.DataFrame, nn: pd.DataFrame) -> dict:
    """Analysis C: BioEmu-1's raw SAXS error against the closeness of each construct to the training data.
    Closeness: the identity of the nearest training sequence (local alignment covering at least half of the
    construct), in four bins and as a continuous variable (Spearman), and whether a CATH superfamily of the
    construct occurs in the training data (any component, or the MD components), present against absent
    (Mann-Whitney; constructs without a CATH superfamily form a third group, reported but not tested).
    Pooled and within each disorder class; Holm adjustment within each family of within-class tests; the
    class-stratified van Elteren test for identity >= 50% and for superfamily presence."""
    m = _bioemu_with_training(d, nn)
    out = {"definition": "nearest training identity: highest MMseqs2 identity to any BioEmu-1 training "
                         "sequence (AFDB cluster members, CATH1/CATH2 and MegaSim MD systems, MEGAscale wild "
                         "types, other public MD proteins) over an alignment covering at least half of the "
                         "construct; 0 when no such hit. Superfamily presence: a CATH-Gene3D superfamily "
                         "inside the construct occurs in the training data.",
           "n": int(len(m)), "bins": ID_LABELS}
    seed = SEED + 500
    out["pooled"], out["by_class"] = {}, {}
    for b in ID_LABELS:
        out["pooled"][b] = summarise(m[m["identity_bin"] == b], seed)
        seed += 10
    for cls in analyse.CLASS_ORDER:
        g = m[m["disorder_class"] == cls]
        out["by_class"][cls] = {}
        for b in ID_LABELS:
            out["by_class"][cls][b] = summarise(g[g["identity_bin"] == b], seed)
            seed += 10
    # identity: Spearman, pooled and within class
    out["spearman_pooled"] = {k: _spearman(m["nearest_identity"], m[k]) for k in CLOSE_METRICS}
    raw = {}
    out["spearman_by_class"] = {}
    for cls in analyse.CLASS_ORDER:
        g = m[m["disorder_class"] == cls]
        out["spearman_by_class"][cls] = {k: _spearman(g["nearest_identity"], g[k]) for k in CLOSE_METRICS}
        for k in CLOSE_METRICS:
            if out["spearman_by_class"][cls][k]["p"] is not None:
                raw[f"{cls}|{k}"] = out["spearman_by_class"][cls][k]["p"]
    for key, p in holm(raw).items():
        cls, k = key.split("|")
        out["spearman_by_class"][cls][k]["holm_p"] = p
    out["holm_family_identity"] = f"the {len(raw)} within-class Spearman tests (three classes, four measures)"
    k50 = m["nearest_identity"] >= 0.5
    out["stratified_identity_ge50"] = {
        "chi2_raw": van_elteren(np.log(m["chi2_raw"].to_numpy(float)), k50.to_numpy(),
                                m["disorder_class"].to_numpy()),
        "n_ge50": int(k50.sum())}
    # superfamily presence
    for grp in ("sf_group", "sf_group_md"):
        rec = {"pooled": {}, "by_class": {}}
        for gname in ("present", "absent", "none"):
            rec["pooled"][gname] = summarise(m[m[grp] == gname], seed)
            seed += 10
        a, b = m[m[grp] == "present"], m[m[grp] == "absent"]
        rec["pooled_tests"] = {f"{k}_mannwhitney_p": _mw(a[k], b[k]) for k in MW_METRICS}
        raw = {}
        for cls in analyse.CLASS_ORDER:
            g = m[m["disorder_class"] == cls]
            r = {gname: summarise(g[g[grp] == gname], seed + i)
                 for i, gname in enumerate(("present", "absent", "none"))}
            seed += 10
            a, b = g[g[grp] == "present"], g[g[grp] == "absent"]
            r["tests"] = {f"{k}_mannwhitney_p": _mw(a[k], b[k]) for k in MW_METRICS}
            for k, p in r["tests"].items():
                if not np.isnan(p):
                    raw[f"{cls}|{k}"] = p
            rec["by_class"][cls] = r
        for key, p in holm(raw).items():
            cls, k = key.split("|")
            rec["by_class"][cls]["tests"][k.replace("_mannwhitney_p", "_holm_p")] = p
        rec["holm_family"] = f"the {len(raw)} within-class Mann-Whitney tests with both groups present"
        kk = m[grp].isin(["present", "absent"])
        rec["stratified_by_class"] = {"chi2_raw": van_elteren(np.log(m.loc[kk, "chi2_raw"].to_numpy(float)),
                                                              (m.loc[kk, grp] == "present").to_numpy(),
                                                              m.loc[kk, "disorder_class"].to_numpy())}
        out[grp] = rec
    out["class_by_bin"] = {b: {c: int(((m["identity_bin"] == b) & (m["disorder_class"] == c)).sum())
                               for c in analyse.CLASS_ORDER} for b in ID_LABELS}
    return out


ALPHA_BETA_ENZYME = ("3.40.50.", "3.20.20.")


def category_evidence(d: pd.DataFrame, mp: pd.DataFrame, nn: pd.DataFrame, rk: pd.DataFrame) -> dict:
    """SAXS-archive side of the candidate categories for data collection: for each category, the clean
    profiles in it and BioEmu-1's raw error on them (raw-fit share, medians of raw chi2, NRMSD and Rg ratio),
    with the same numbers for the rest of the set. The training-data side is in training_summary.json."""
    m = _bioemu_with_training(d, nn)
    mpi = mp.set_index("label")
    t = nn.set_index("label")
    lab = m["label"]
    sfl = lab.map(mpi["gene3d_superfamilies"]).fillna("").str.split(";").map(lambda x: [y for y in x if y])
    rank = dict(zip(rk["sf"], rk["rank"], strict=True))
    grp = structure_group(mpi)
    cats = {
        "disordered": m["disorder_class"] == "disordered",
        "partly disordered": m["disorder_class"] == "partly disordered",
        "mostly disordered (metapredict)": lab.map(t["mp_dis_cat"]) == "mostly disordered",
        "long disordered region, under half": lab.map(t["mp_dis_cat"]) == "long disordered region",
        "construct over 500 residues": m["length"] > 500,
        "two or more CATH superfamilies": sfl.map(lambda x: len(set(x)) >= 2),
        "alpha/beta enzyme folds": sfl.map(lambda x: any(y.startswith(ALPHA_BETA_ENZYME) for y in x)),
        "mainly-beta domains": sfl.map(lambda x: any(y.startswith("2.") for y in x)),
        "bacteria other than E. coli": lab.map(mpi["org_group"]) == "Other bacteria",
        "human proteins": lab.map(mpi["org_group"]) == "Human",
        "superfamily beyond the top 200": sfl.map(lambda x: any(rank.get(y, 10 ** 6) > 200 for y in x)),
        "no experimental structure": lab.map(grp) == "none",
    }
    out = {"n": int(len(m)), "categories": {}}
    seed = SEED + 900
    for name, k in cats.items():
        k = k.fillna(False).astype(bool).to_numpy()
        out["categories"][name] = {"in": summarise(m[k], seed), "rest": summarise(m[~k], seed + 5),
                                   "chi2_raw_mannwhitney_p": _mw(m.loc[k, "chi2_raw"], m.loc[~k, "chi2_raw"]),
                                   "rg_ratio_mannwhitney_p": _mw(m.loc[k, "rg_ratio"], m.loc[~k, "rg_ratio"])}
        seed += 10
    out["max_construct_length"] = int(m["length"].max())
    return out


def quality_checks_by_class(ent: pd.DataFrame) -> dict:
    """Share of all archive entries that pass the quality checks (analyse.clean_mask), by disorder class, with
    Wilson intervals, the count stopped by each check, and a two-sided Fisher exact test of disordered against
    folded. An entry can raise more than one flag, so per-check counts can sum to more than the failures."""
    ok = analyse.clean_mask(ent)
    checks = {"aggregation upturn": ent["flag_aggregation"].astype(bool),
              "negative intensities": ent["flag_negative_I"].astype(bool),
              "unusable errors": ent["flag_bad_errors"].astype(bool),
              "mass above a monomer": ent["flag_oligomer"].astype(bool),
              "no valid Guinier region": ~ent["guinier_valid"].astype(bool)}
    out = {"n": int(len(ent)), "rule": "no data-quality flag raised and a valid Guinier region",
           "by_class": {}}
    for c in analyse.CLASS_ORDER:
        k = (ent["disorder_class"] == c).to_numpy()
        n, kp = int(k.sum()), int(ok[k].sum())
        out["by_class"][c] = {"n": n, "pass": kp, "share": kp / n, "ci95": wilson(kp, n),
                              "failed_by_check": {nm: int(v[k].sum()) for nm, v in checks.items()}}
    d_, f_ = out["by_class"]["disordered"], out["by_class"]["folded"]
    table = [[d_["pass"], d_["n"] - d_["pass"]], [f_["pass"], f_["n"] - f_["pass"]]]
    out["fisher_disordered_vs_folded"] = {"table_pass_fail": table,
                                          "odds_ratio": float(stats.fisher_exact(table)[0]),
                                          "p_two_sided": float(stats.fisher_exact(table)[1])}
    return out


# ----------------------------------------------------------------------------- figures
SAXS_COL = report.PAGE_ACCENT       # the SAXS set
REF_COL = report.NEUTRAL_BAR        # the experimental reference (CATH S35, PDB clusters)
SP_COL = "#C3C9CF"                  # Swiss-Prot


def fig_coverage_space(S: dict, rc: dict) -> None:
    """A: share of CATH superfamilies present in the SAXS set by rank of their S35 cluster count.
    B: CATH class composition of the SAXS set against CATH S35 clusters."""
    fig, axes = report.figure(720, 300, ncols=2, gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    b = rc["bins"]
    x = np.arange(len(b))
    sh = np.array([r["share"] for r in b])
    lo = sh - np.array([r["ci95"][0] for r in b])
    hi = np.array([r["ci95"][1] for r in b]) - sh
    ax.bar(x, sh, width=0.68, color=SAXS_COL)
    ax.errorbar(x, sh, yerr=[lo, hi], fmt="none", ecolor=report.INK, lw=0.9, capsize=2.5)
    for i, r in enumerate(b):
        ax.text(i, r["ci95"][1] + 0.025, f"{r['k']}/{r['n']}", ha="center", fontsize=9, color=report.MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['first']}–\n{r['last']}" for r in b], fontsize=9.5)
    ax.set_xlabel("CATH superfamilies, ranked by number of S35 clusters")
    ax.set_ylabel("share present in the SAXS set")
    ax.set_ylim(0, 1.0)
    ax.set_title("A · common superfamilies are covered, rare ones are not", loc="left")
    ax = axes[1]
    cls = pd.DataFrame(S["composition"]["class"])
    cls = cls[cls["class"].isin([1, 2, 3])]
    y = np.arange(len(cls))
    h = 0.36
    ax.barh(y - h / 2 - 0.02, cls["ours_share"], h, color=SAXS_COL, label="SAXS set")
    ax.errorbar(cls["ours_share"], y - h / 2 - 0.02, xerr=[cls["ours_share"] - cls["ours_lo"],
                                                         cls["ours_hi"] - cls["ours_share"]],
                fmt="none", ecolor=report.INK, lw=0.9, capsize=2.5)
    ax.barh(y + h / 2 + 0.02, cls["cath_share"], h, color=REF_COL,
            label=f"CATH {S['cath_reference']['version']} S35")
    for yy, (_, r) in zip(y, cls.iterrows(), strict=True):
        ax.text(max(r["ours_hi"], r["cath_share"]) + 0.015, yy, f"×{r['enrichment']:.2f}", va="center",
                fontsize=9.5, color=report.INK)
    ax.set_yticks(y)
    ax.set_yticklabels(["mainly α", "mainly β", "α/β"])
    ax.invert_yaxis()
    ax.set_xlim(0, 0.86)
    ax.set_xlabel("share of protein–superfamily units")
    ax.legend(frameon=False, loc="upper right", fontsize=9.5)
    ax.set_title("B · CATH class (× = SAXS share / CATH share)", loc="left")
    fig.tight_layout(w_pad=2.0)
    report.save(fig, "fig_coverage_space")


def fig_coverage_sources(S: dict) -> None:
    """A: organism of the SAXS set against PDB 30% clusters and Swiss-Prot. B: CATH-Gene3D domains
    in the construct, in the full-length parent protein and in Swiss-Prot proteins."""
    fig, axes = report.figure(720, 300, ncols=2, gridspec_kw={"width_ratios": [1.1, 1]})
    ax = axes[0]
    org = pd.DataFrame(S["organism"]["groups"]).set_index("group")
    order = ["Human", "Mouse", "E. coli", "Yeast", "Other bacteria", "Viruses", "Plants", "Other"]
    org = org.loc[order]
    y = np.arange(len(org))
    h = 0.26
    ax.barh(y - h - 0.02, org["ours_share"], h, color=SAXS_COL, label="SAXS set")
    ax.errorbar(org["ours_share"], y - h - 0.02,
                xerr=[org["ours_share"] - org["ours_lo"], org["ours_hi"] - org["ours_share"]],
                fmt="none", ecolor=report.INK, lw=0.8, capsize=2)
    ax.barh(y, org["pdb_nr30_share"], h, color=REF_COL, label="PDB, 30% clusters")
    ax.barh(y + h + 0.02, org["swissprot_share"], h, color=SP_COL, label="Swiss-Prot")
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.6)
    ax.set_xlabel("share of proteins")
    ax.legend(frameon=False, loc="lower right", fontsize=9.5)
    ax.set_title("A · source organism", loc="left")
    ax = axes[1]
    dp = S["domains_per_protein"]
    cats = ["0", "1", "2", "3", "4+"]
    x = np.arange(len(cats))
    w = 0.27
    for off, key, col, lab in ((-w - 0.02, "ours_construct_gene3d", SAXS_COL, "SAXS construct"),
                               (0.0, "ours_parent_fulllength_gene3d", "#9DB4CB", "full-length protein"),
                               (w + 0.02, "swissprot_fulllength_gene3d", SP_COL, "Swiss-Prot")):
        ax.bar(x + off, [dp[key][c] for c in cats], w, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_xlabel("CATH-Gene3D domains")
    ax.set_ylabel("share of proteins")
    ax.set_ylim(0, 0.46)
    ax.legend(frameon=False, loc="upper right", fontsize=9.5)
    ax.set_title("B · domains per protein", loc="left")
    fig.tight_layout(w_pad=2.0)
    report.save(fig, "fig_coverage_sources")


def fig_coverage_architecture(S: dict) -> None:
    """Enrichment of CATH architectures in the SAXS set over CATH S35, with bootstrap intervals."""
    ar = pd.DataFrame(S["composition"]["architecture"]).sort_values("enrichment")
    fig, ax = report.figure(640, 420)
    y = np.arange(len(ar))
    e = ar["enrichment"].clip(lower=0.1)
    ax.axvline(1, color=report.MUTED, lw=0.8, ls="--")
    ax.errorbar(e, y, xerr=[e - ar["enr_lo"].clip(lower=0.1), ar["enr_hi"] - e], fmt="o", color=SAXS_COL,
                ecolor=SAXS_COL, ms=5, capsize=2, lw=1.1)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 8)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4, 8])
    ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "8"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['arch']} {r['name']} ({r['ours_units']})" for _, r in ar.iterrows()],
                       fontsize=9.5)
    ax.set_xlabel("SAXS share / CATH S35 share (log scale, 95% bootstrap interval)")
    ax.set_title("CATH architectures (CATH share ≥ 1% or ≥ 3 units in the SAXS set)", loc="left")
    fig.tight_layout()
    report.save(fig, "fig_coverage_architecture")


GROUP_ALPHA = {"pdb50": 0.75, "partial": 0.45, "none": 0.75}


def fig_structure_split(d: pd.DataFrame, mp: pd.DataFrame, SP: dict) -> None:
    """Raw chi2 of BioEmu-1 by structure group within each disorder class: points, the group median and
    its 95% interval; filled circles PDB >= 50%, light circles partial, open circles none."""
    m = d[(d["model"] == C.MODEL) & d["clean"]].copy()
    m["structure"] = m["label"].map(structure_group(mp.set_index("label")))
    m = m.dropna(subset=["structure"])
    fig, ax = report.figure(720, 330)
    ticks, labels = [], []
    pos = 0.0
    for ci_, cls in enumerate(analyse.CLASS_ORDER):
        col = report.CLS_COL[cls]
        for gi, gname in enumerate(GROUPS):
            sel = (m["disorder_class"] == cls) & (m["structure"] == gname)
            v = np.log10(m.loc[sel, "chi2_raw"].to_numpy())
            x = pos + np.random.default_rng(100 * ci_ + gi).uniform(-0.2, 0.2, len(v))
            if gname == "none":
                ax.scatter(x, v, s=13, facecolors="none", edgecolors=col, linewidths=0.8, alpha=0.8)
            else:
                ax.scatter(x, v, s=13, color=col, alpha=GROUP_ALPHA[gname], lw=0)
            st = SP["by_class"][cls][gname]
            if st["n"]:
                med = np.log10(st["chi2_raw_median"])
                ax.plot([pos - 0.28, pos + 0.28], [med, med], color=report.INK, lw=2)
                if st["n"] >= 5:
                    ax.plot([pos + 0.33] * 2, np.log10(st["chi2_raw_ci95"]), color=report.INK, lw=1.4)
            ticks.append(pos)
            labels.append(f"{GROUP_LABEL[gname]}\nn = {st['n']}")
            pos += 1
        t = SP["by_class"][cls]["tests"]
        ax.text(pos - 2, 3.55, f"{cls}", ha="center", fontsize=10, color=report.INK, fontweight="bold")
        ax.text(pos - 2, 3.3, f"Mann–Whitney p = {t['chi2_raw_mannwhitney_p']:.2f}", ha="center",
                fontsize=9.5, color=report.MUTED)
        pos += 0.7
    ax.axhline(np.log10(RAW_FIT), color=report.MUTED, lw=0.8, ls="--")
    ax.text(pos - 0.75, np.log10(RAW_FIT) - 0.05, "χ² = 2", fontsize=9.5, color=report.MUTED, ha="right",
            va="top")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_xlim(-0.6, pos - 0.6)
    ax.set_ylim(-0.6, 3.75)
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["1", "10", "100", "1000"])
    ax.set_ylabel("raw reduced χ² (log scale)")
    fig.tight_layout()
    report.save(fig, "fig_structure_split")


# training-data figures: one colour per set (validated for colour-vision deficiency, all pairs); the CATH and
# PDB references stay in the light neutral of the earlier coverage figures
SET_COL = {"SAXS": SAXS_COL, "AFDB": "#139C87", "CATH MD": "#C77D12", "CATH S35": SP_COL, "PDB 30%": SP_COL}
SET_LABEL = {"SAXS": "SAXS set", "AFDB": "AFDB pre-training clusters", "CATH MD": "CATH MD domains",
             "CATH S35": "CATH S35 (reference)", "PDB 30%": "PDB 30% clusters (reference)"}


def _grouped_barh(ax, cats, sets, value, h_total=0.8):
    """Horizontal grouped bars: one group per category, one bar per set (top to bottom in set order)."""
    n = len(sets)
    h = h_total / n
    y = np.arange(len(cats))
    for j, st in enumerate(sets):
        off = -h_total / 2 + h * (j + 0.5)
        ax.barh(y + off, [value(st, c) for c in cats], h * 0.86, color=SET_COL[st], label=SET_LABEL[st])
    ax.set_yticks(y)
    ax.invert_yaxis()
    return y


def fig_training_composition(T: dict) -> None:
    """Training data against the SAXS set and the references: A CATH class, B chain length, C predicted
    disorder, D organism."""
    fig, axes = report.figure(760, 600, nrows=2, ncols=2)
    sets4 = ["SAXS", "AFDB", "CATH MD", "CATH S35"]
    ax = axes[0, 0]
    cls = ["Mainly alpha", "Mainly beta", "Alpha beta"]
    _grouped_barh(ax, cls, sets4, lambda st, c: T["composition"][st]["class"][c]["share"])
    ax.set_yticklabels(["mainly α", "mainly β", "α/β"])
    ax.set_xlim(0, 0.62)
    ax.set_xlabel("share of protein–superfamily units")
    ax.set_title("A · CATH class", loc="left")
    ax = axes[0, 1]
    for j, st in enumerate(sets4):
        L = T["length"][st]
        ax.plot([L["p10"], L["p90"]], [j, j], color=SET_COL[st], lw=1.6)
        ax.plot([L["p25"], L["p75"]], [j, j], color=SET_COL[st], lw=7, solid_capstyle="butt")
        ax.plot([L["median"]] * 2, [j - 0.26, j + 0.26], color=report.INK, lw=1.6)
        ax.text(L["p90"] * 1.06, j, f"{L['median']:.0f}", va="center", fontsize=9.5, color=report.MUTED)
    ax.axvline(500, color=report.MUTED, lw=0.8, ls="--")
    ax.set_xscale("log")
    ax.set_xlim(25, 2500)
    ax.set_xticks([30, 50, 100, 200, 500, 1000, 2000])
    ax.set_xticklabels(["30", "50", "100", "200", "500", "1000", "2000"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(range(len(sets4)))
    ax.set_yticklabels(["SAXS constructs", "AFDB inputs", "CATH MD", "CATH S35"])
    ax.invert_yaxis()
    ax.set_xlabel("chain length (residues, log scale)")
    ax.set_title("B · length: median, 25–75% and 10–90%", loc="left")
    ax = axes[1, 0]
    dis = ["mostly disordered", "long disordered region"]
    _grouped_barh(ax, dis, sets4, lambda st, c: T["disorder"][st][c]["share"])
    ax.set_yticklabels(["≥ 50% of residues\ndisordered", "a disordered run\nof ≥ 30, under 50%"])
    ax.set_xlabel("share of chains (metapredict)")
    ax.set_xlim(0, max(T["disorder"][st][c]["share"] for st in sets4 for c in dis) * 1.15)
    ax.set_title("C · predicted disorder", loc="left")
    ax = axes[1, 1]
    orgs = ["Human", "Mouse", "E. coli", "Yeast", "Other bacteria", "Viruses", "Plants", "Other"]
    sets_o = ["SAXS", "AFDB", "CATH MD", "PDB 30%"]
    _grouped_barh(ax, orgs, sets_o, lambda st, c: T["organism"][st][c]["share"], h_total=0.84)
    ax.set_yticklabels(orgs, fontsize=9)
    ax.set_xlim(0, 0.95)
    v = T["organism"]["AFDB"]["Other bacteria"]["share"]       # the AFDB bar, second of four in its group
    ax.text(v + 0.01, orgs.index("Other bacteria") - 0.105, f"{v:.0%}", va="center", fontsize=9,
            color=report.MUTED)
    ax.set_xlabel("share of chains")
    ax.set_title("D · organism", loc="left")
    handles = [matplotlib.patches.Patch(color=SET_COL[st], label=SET_LABEL[st])
               for st in ["SAXS", "AFDB", "CATH MD", "CATH S35"]]
    handles[-1] = matplotlib.patches.Patch(color=SP_COL,
                                           label="reference: CATH S35 (A–C), PDB 30% clusters (D)")
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=9.5)
    fig.tight_layout(rect=(0, 0.075, 1, 1), h_pad=1.6, w_pad=2.0)
    report.save(fig, "fig_training_composition")


def fig_training_overlap(nn: pd.DataFrame) -> None:
    """A: share of the clean constructs whose nearest training sequence reaches each identity, for all
    training data, the AFDB clusters and the MD sets. B: the constructs by identity bin and disorder class."""
    c = nn[nn["clean"].astype(str) == "True"]
    fig, axes = report.figure(760, 300, ncols=2, gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    x = np.linspace(0, 1, 201)
    for col, lab, colour, lw in (("nearest_identity", "any training data", SAXS_COL, 2.2),
                                 ("id_afdb_member", "AFDB clusters", SET_COL["AFDB"], 1.6),
                                 ("nearest_identity_md", "MD sets (CATH, MegaSim)", SET_COL["CATH MD"], 1.6)):
        v = c[col].to_numpy(float)
        ax.plot(x, [(v >= t - 1e-9).mean() for t in x], color=colour, lw=lw, label=lab,
                alpha=report.LINE_ALPHA)
    for t in (0.3, 0.5, 0.9):
        ax.axvline(t, color=report.MUTED, lw=0.7, ls="--")
    ax.set_xlim(0.2, 1.0)
    ax.set_ylim(0, 0.7)
    ax.set_xticks([0.2, 0.3, 0.5, 0.9, 1.0])
    ax.set_xticklabels(["20%", "30%", "50%", "90%", "100%"])
    ax.set_xlabel("identity of the nearest training sequence")
    ax.set_ylabel("share of constructs at or above")
    ax.legend(frameon=False, loc="upper right", fontsize=9.5)
    ax.set_title("A · closeness to the training data", loc="left")
    ax = axes[1]
    b = identity_bin(c["nearest_identity"])
    left = np.zeros(len(ID_LABELS))
    for cls in analyse.CLASS_ORDER:
        v = np.array([((b == lab) & (c["disorder_class"] == cls)).sum() for lab in ID_LABELS])
        ax.barh(np.arange(len(ID_LABELS)), v, left=left, color=report.CLS_COL[cls], label=cls, height=0.62,
                edgecolor="white", lw=1)
        left += v
    for i, n in enumerate(left):
        ax.text(n + 3, i, f"{int(n)}", va="center", fontsize=9.5, color=report.MUTED)
    ax.set_yticks(range(len(ID_LABELS)))
    ax.set_yticklabels(["< 30%", "30–50%", "50–90%", "≥ 90%"])
    ax.invert_yaxis()
    ax.set_xlim(0, left.max() * 1.18)
    ax.set_xlabel("constructs")
    ax.set_ylabel("nearest training identity")
    ax.legend(frameon=False, loc="lower right", fontsize=9.5)
    ax.set_title("B · by disorder class", loc="left")
    fig.tight_layout(w_pad=2.0)
    report.save(fig, "fig_training_overlap")


def fig_closeness(d: pd.DataFrame, nn: pd.DataFrame) -> None:
    """Raw chi2 (left) and Rg ratio (right) of BioEmu-1 by nearest-training identity bin, within each
    disorder class; bars are group medians (groups of three or more)."""
    m = _bioemu_with_training(d, nn)
    fig, axes = report.figure(760, 560, nrows=2)
    for ax, col, ylab, ref in ((axes[0], "chi2_raw", "raw reduced χ²\n(log scale)", RAW_FIT),
                               (axes[1], "rg_ratio", "ensemble Rg /\nmeasured Rg", 1.0)):
        pos, ticks, labels = 0.0, [], []
        for ci_, cls in enumerate(analyse.CLASS_ORDER):
            col_c = report.CLS_COL[cls]
            for bi, b in enumerate(ID_LABELS):
                v = m.loc[(m["disorder_class"] == cls) & (m["identity_bin"] == b), col].dropna().to_numpy()
                if col == "chi2_raw":
                    v = np.log10(v)
                x = pos + np.random.default_rng(10 * ci_ + bi).uniform(-0.18, 0.18, len(v))
                ax.scatter(x, v, s=9, color=col_c, alpha=0.55, lw=0)
                if len(v) >= 3:
                    ax.plot([pos - 0.3, pos + 0.3], [np.median(v)] * 2, color=report.INK, lw=1.8)
                ticks.append(pos)
                labels.append(b.replace(">=", "≥ ").replace("<", "< ").replace("-", "–"))
                pos += 1
            ax.text(pos - 2.5, 1.0, cls, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=9.5, color=report.INK)
            pos += 0.8
        ax.axhline(np.log10(ref) if col == "chi2_raw" else ref, color=report.MUTED, lw=0.8, ls="--")
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_xlim(-0.6, pos - 0.6)
        ax.set_ylabel(ylab)
        if col == "rg_ratio":
            ax.set_xlabel("identity of the nearest training sequence (%)")
        if col == "chi2_raw":
            ax.set_yticks([0, 1, 2, 3])
            ax.set_yticklabels(["1", "10", "100", "1000"])
        else:
            ax.set_ylim(0.5, 2.0)
    fig.tight_layout(w_pad=2.0)
    report.save(fig, "fig_closeness")


def main() -> None:
    report.style()
    report.OUT.mkdir(parents=True, exist_ok=True)
    mp = pd.read_csv(COV / "per_protein_mapping.csv")
    S = json.load(open(COV / "summary.json"))
    rk = pd.read_csv(COV / "cath_superfamily_rank.csv.gz")
    d = pd.read_csv(C.RESULTS / "resolvability.csv")
    SP = structure_split(d, mp)
    json.dump(SP, open(COV / "structure_split.json", "w"), indent=1)
    rc = rank_coverage(rk)
    json.dump(rc, open(COV / "rank_coverage.json", "w"), indent=1)
    fig_coverage_space(S, rc)
    fig_coverage_sources(S)
    fig_coverage_architecture(S)
    fig_structure_split(d, mp, SP)
    nn = pd.read_csv(COV / "training_nearest.csv")
    T = json.load(open(COV / "training_summary.json"))
    TC = training_closeness(d, nn)
    json.dump(TC, open(COV / "training_closeness.json", "w"), indent=1)
    CE = category_evidence(d, mp, nn, rk)
    CE["quality_checks_by_class"] = quality_checks_by_class(pd.read_csv(C.RESULTS / "entries.csv"))
    json.dump(CE, open(COV / "category_evidence.json", "w"), indent=1)
    fig_training_composition(T)
    fig_training_overlap(nn)
    fig_closeness(d, nn)
    print(f"  structure split and rank coverage -> {COV}; figures -> {report.OUT}", flush=True)


if __name__ == "__main__":
    main()
