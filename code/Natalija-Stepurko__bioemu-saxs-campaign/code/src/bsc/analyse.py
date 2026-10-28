"""Stage `analyse`: the error map by protein class, the four pre-specified expectations, the
reweighting kind of every entry, model comparisons, and the ranked list of candidate systems for
new measurements; with the robustness statistics added after the first run (DESIGN §7).

Output: results/analysis.json (every number the page reports), results/by_class.csv,
results/resolvability.csv, results/candidates.csv, results/clusters.csv (sequence groups).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import stats

from bsc import config as C
from bsc import robust

CLASS_ORDER = [c[0] for c in C.DISORDER_CLASSES]
# the four kinds along the reweighting path (renamed after the first run and again on 2026-10-09,
# DESIGN §7; the logic is unchanged: raw chi2 <= 2; else chi2 <= 2 reached with phi >= threshold;
# else reached; else never)
KIND_RAW, KIND_MODEST = "raw fit", "modestly reweightable"
KIND_STRONG, KIND_NOTFIT = "strongly reweightable", "target fit not reached"
KINDS = [KIND_RAW, KIND_MODEST, KIND_STRONG, KIND_NOTFIT]
PHI_THRESHOLD = 0.5
PHI_SENSITIVITY = (0.3, 0.5, 0.7)
N_BOOT = 2000
BOOT_SEED = 20261009


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    ent = pd.read_csv(C.RESULTS / "entries.csv")
    sc = pd.read_csv(C.RESULTS / "scores.csv")
    df = sc.merge(ent.drop(columns=["sequence"]), on="label", how="left")
    df["log_chi2_raw"] = np.log10(df["chi2_raw"])
    df["rg_ratio"] = df["rg_model"] / df["rg_exp"]
    df["clean"] = clean_mask(df)
    return ent, df


QC_FLAGS = ["flag_aggregation", "flag_negative_I", "flag_bad_errors", "flag_oligomer"]


def clean_mask(df: pd.DataFrame) -> pd.Series:
    """An entry passes the quality checks when no data-quality flag is raised and the Guinier fit is valid."""
    flagged = df[QC_FLAGS].astype(bool).any(axis=1)
    return ~flagged & df["guinier_valid"].astype(bool)


def kind(row: pd.Series, phi_threshold: float = PHI_THRESHOLD) -> str:
    """One label per entry from the raw fit and the reweighting path (thresholds fixed in DESIGN §4).

    raw fit                   raw chi2 <= 2: the prior ensemble already describes the data
    modestly reweightable     chi2 <= 2 reached while keeping phi >= threshold (0.5)
    strongly reweightable     chi2 <= 2 reached only with phi < threshold
    target fit not reached    chi2 <= 2 never reached along the path
    """
    if row["chi2_raw"] <= 2:
        return KIND_RAW
    if np.isfinite(row["phi_at_chi2_2"]) and row["phi_at_chi2_2"] >= phi_threshold:
        return KIND_MODEST
    if np.isfinite(row["phi_at_chi2_2"]):
        return KIND_STRONG
    return KIND_NOTFIT


resolvability = kind      # name used by the first run


def choose_examples(d: pd.DataFrame) -> dict[str, str]:
    """The worked example of each class: the clean entry of the model under study whose raw chi2 is
    closest to the class median (an even count leaves two entries equidistant; the shorter chain is
    taken, then the alphabetically first label). A joint rule on raw chi2 and size ratio was tried and
    not kept (DESIGN §7)."""
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    out = {}
    for cls in CLASS_ORDER:
        g = m[m["disorder_class"] == cls].copy()
        g["dev"] = (g["chi2_raw"] - g["chi2_raw"].median()).abs().round(9)
        out[cls] = str(g.sort_values(["dev", "length", "label"]).iloc[0]["label"])
    return out


def choose_representative_path(d: pd.DataFrame) -> str | None:
    """The one reweighting path drawn on the page: the strongly reweightable clean entry whose phi at
    chi2 = 2 is closest to the median of that kind (ties: shorter chain, then label)."""
    m = d[(d["model"] == C.MODEL) & d["clean"] & (d["kind"] == KIND_STRONG)].copy()
    if m.empty:
        return None
    m["dev"] = (m["phi_at_chi2_2"] - m["phi_at_chi2_2"].median()).abs().round(9)
    return str(m.sort_values(["dev", "length", "label"]).iloc[0]["label"])


def ci(values, stat=np.median, seed_offset: int = 0) -> list[float]:
    lo, hi = robust.bootstrap_ci(np.asarray(values, float), stat, N_BOOT, BOOT_SEED + seed_offset)
    return [lo, hi]


def expectations(d: pd.DataFrame) -> dict:
    """E1–E4 of DESIGN §5 on the model under study, clean entries only."""
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    out = {}
    f, dis = m[m["disorder_class"] == "folded"], m[m["disorder_class"] == "disordered"]
    u = stats.mannwhitneyu(dis["log_chi2_raw"], f["log_chi2_raw"], alternative="greater")
    out["E1"] = {"statement": "raw chi2 worse for disordered than folded",
                 "median_folded": float(f["chi2_raw"].median()),
                 "median_disordered": float(dis["chi2_raw"].median()),
                 "n_folded": len(f), "n_disordered": len(dis), "p_one_sided": float(u.pvalue),
                 "met": bool(dis["chi2_raw"].median() > f["chi2_raw"].median() and u.pvalue < 0.05)}
    rho = stats.spearmanr(f["length"], f["log_chi2_raw"])
    small, large = f[f["length"] <= 100], f[f["length"] > 100]
    out["E2"] = {"statement": "among folded proteins raw chi2 rises with length",
                 "spearman_rho": float(rho.statistic), "p": float(rho.pvalue), "n": len(f),
                 "share_fit_le_100": float((small["chi2_raw"] <= 2).mean()), "n_le_100": len(small),
                 "share_fit_gt_100": float((large["chi2_raw"] <= 2).mean()), "n_gt_100": len(large),
                 "median_chi2_le_100": float(small["chi2_raw"].median()),
                 "median_chi2_gt_100": float(large["chi2_raw"].median()),
                 "met": bool(rho.statistic > 0 and rho.pvalue < 0.05)}
    reached = m["phi_at_chi2_1"].notna()
    out["E3"] = {"statement": "most entries reach chi2 ~1 with phi above 0.3; a minority need phi below 0.1",
                 "share_reaching_chi2_1": float(reached.mean()),
                 "share_phi_above_0.3": float((m["phi_at_chi2_1"] > 0.3).mean()),
                 "share_phi_below_0.1": float((m["phi_at_chi2_1"] < 0.1).mean()),
                 "met": bool((m["phi_at_chi2_1"] > 0.3).mean() > 0.5)}
    ratio = dis["rg_ratio"].dropna()
    logr = np.log(ratio) if len(ratio) > 5 else np.array([0, 0, 0, 0, 0, 0])
    w = stats.wilcoxon(logr, alternative="less")
    w2 = stats.wilcoxon(logr, alternative="two-sided")
    out["E4"] = {"statement": "ensemble Rg below experimental Rg for disordered proteins",
                 "median_rg_ratio_disordered": float(ratio.median()), "n": int(len(ratio)),
                 "share_below_1": float((ratio < 1).mean()), "p_one_sided": float(w.pvalue),
                 "p_two_sided": float(w2.pvalue), "ci95_median_rg_ratio": ci(ratio, seed_offset=4),
                 "met": bool(ratio.median() < 1 and w.pvalue < 0.05)}
    return out


def by_class(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, cls), g in d[d["clean"]].groupby(["model", "disorder_class"]):
        rows.append({"model": model, "disorder_class": cls, "n": len(g),
                     "chi2_raw_median": g["chi2_raw"].median(),
                     "chi2_raw_q25": g["chi2_raw"].quantile(0.25),
                     "chi2_raw_q75": g["chi2_raw"].quantile(0.75),
                     "share_raw_fit": (g["chi2_raw"] <= 2).mean(),
                     "chi2_best_median": g["chi2_best"].median(),
                     "phi_at_chi2_1_median": g["phi_at_chi2_1"].median(),
                     "share_not_fit": (g["kind"] == KIND_NOTFIT).mean(),
                     "rg_ratio_median": g["rg_ratio"].median(),
                     "nrmsd_raw_median": g["nrmsd_raw"].median() if "nrmsd_raw" in g else np.nan,
                     "nrmsd_best_median": g["nrmsd_best"].median() if "nrmsd_best" in g else np.nan})
    out = pd.DataFrame(rows)
    out["disorder_class"] = pd.Categorical(out["disorder_class"], CLASS_ORDER, ordered=True)
    return out.sort_values(["model", "disorder_class"])


def model_comparison(d: pd.DataFrame) -> list[dict]:
    """Paired comparison of each comparator with the model under study on the entries both have,
    with bootstrap intervals over entries for the median log10 chi2 ratio and the share better."""
    base = d[(d["model"] == C.MODEL) & d["clean"]].set_index("label")
    out = []
    for k, m in enumerate(sorted(d["model"].unique())):
        if m == C.MODEL:
            continue
        o = d[(d["model"] == m) & d["clean"]].set_index("label")
        common = base.index.intersection(o.index)
        if len(common) < 10:
            continue
        diff = (np.log10(o.loc[common, "chi2_raw"]) - np.log10(base.loc[common, "chi2_raw"])).to_numpy()
        w = stats.wilcoxon(diff)
        rec = {"model": m, "n": int(len(common)),
               "median_log10_chi2_ratio_vs_bioemu": float(np.median(diff)),
               "ci95_median_log10_chi2_ratio": ci(diff, seed_offset=100 + k),
               "share_better_than_bioemu": float((diff < 0).mean()),
               "ci95_share_better": ci((diff < 0).astype(float), robust.share, seed_offset=200 + k),
               "p_two_sided": float(w.pvalue),
               "median_chi2_raw": float(o.loc[common, "chi2_raw"].median()),
               "bioemu_median_chi2_raw_same_entries": float(base.loc[common, "chi2_raw"].median())}
        if "nrmsd_raw" in o:
            rec["median_nrmsd_raw"] = float(o.loc[common, "nrmsd_raw"].median())
            rec["bioemu_median_nrmsd_raw_same_entries"] = float(base.loc[common, "nrmsd_raw"].median())
            nd = (o.loc[common, "nrmsd_raw"] - base.loc[common, "nrmsd_raw"]).to_numpy()
            rec["share_better_nrmsd"] = float(np.nanmean(nd < 0))
        out.append(rec)
    return out


def candidates(d: pd.DataFrame) -> pd.DataFrame:
    """Rank entries for new measurement: the model is wrong (large raw chi2 after clean data), the
    data could tell (reweighting reaches a fit), and the system is tractable (length, folded or
    partly disordered). Score = log10 raw chi2 × reweighting weight × tractability weight."""
    m = d[(d["model"] == C.MODEL) & d["clean"]].copy()
    res_w = m["kind"].map({KIND_STRONG: 1.0, KIND_MODEST: 0.6, KIND_NOTFIT: 0.3, KIND_RAW: 0.0})
    tract = np.where(m["length"] <= 350, 1.0, 0.5) * m["disorder_class"].map(
        {"folded": 1.0, "partly disordered": 1.0, "disordered": 0.6}).to_numpy()
    m["priority"] = m["log_chi2_raw"].clip(lower=0) * res_w * tract
    m["rg_direction"] = np.where(m["rg_ratio"] < 0.95, "model too compact",
                                 np.where(m["rg_ratio"] > 1.05, "model too extended", "Rg agrees"))
    m["resolvability"] = m["kind"]
    cols = ["label", "length", "disorder_class", "chi2_raw", "chi2_best", "phi_at_chi2_2",
            "resolvability", "rg_exp", "rg_model", "rg_direction", "mw_ratio", "sasbdb_title", "priority"]
    return m.sort_values("priority", ascending=False)[cols].reset_index(drop=True)


# ---- additions after the first run (DESIGN §7) ---------------------------------------------------
def sequence_clusters(ent: pd.DataFrame) -> pd.DataFrame:
    """Near-duplicate clusters (Jaccard) and the looser containment groups used for cross-validation."""
    seqs = ent["sequence"].tolist()
    dup = robust.cluster_sequences(seqs, robust.NEAR_DUPLICATE, measure=robust.jaccard)
    return pd.DataFrame({"label": ent["label"], "dup_cluster": dup,
                         "group": robust.cluster_sequences(seqs, robust.GROUP_SIMILARITY)})


def _class_ratio_stats(m: pd.DataFrame, col: str = "rg_ratio", seed_offset: int = 0) -> dict:
    out = {}
    for i, cls in enumerate(CLASS_ORDER):
        r = m.loc[m["disorder_class"] == cls, col].dropna()
        if len(r) < 3:
            continue
        rec = {"n": int(len(r)), "median": float(r.median()), "ci95": ci(r, seed_offset=seed_offset + i),
               "share_above_1": float((r > 1).mean())}
        if len(r) > 5:
            rec["p_two_sided_wilcoxon_log"] = float(stats.wilcoxon(np.log(r)).pvalue)
        out[cls] = rec
    return out


def rg_robustness(d: pd.DataFrame, ent: pd.DataFrame, clusters: pd.DataFrame) -> dict:
    """B1: intervals, two-sided test, trimmed sensitivity, length dependence and near-duplicate
    collapse for the Rg ratio finding; B2: coordinate-derived Rg when results/coord_rg.csv exists."""
    m = d[(d["model"] == C.MODEL) & d["clean"]].copy()
    out = {"n_boot": N_BOOT, "seed": BOOT_SEED,
           "curve_rg_ratio_by_class": _class_ratio_stats(m, seed_offset=10)}
    dis = m[m["disorder_class"] == "disordered"].dropna(subset=["rg_ratio"])
    logr = np.log(dis["rg_ratio"])
    keep = logr.abs().rank(ascending=False) > 5
    trimmed = dis[keep]["rg_ratio"]
    out["disordered_drop_5_most_extreme"] = {
        "n": int(len(trimmed)), "median": float(trimmed.median()), "ci95": ci(trimmed, seed_offset=20),
        "p_two_sided_wilcoxon_log": float(stats.wilcoxon(np.log(trimmed)).pvalue),
        "dropped": sorted(dis.loc[~keep, "label"].tolist())}
    rho = stats.spearmanr(dis["length"], dis["rg_ratio"])
    out["disordered_ratio_vs_length"] = {"spearman_rho": float(rho.statistic), "p": float(rho.pvalue),
                                         "n": int(len(dis))}
    pairs = robust.near_duplicate_pairs(ent["label"].tolist(), ent["sequence"].tolist())
    cl = clusters.set_index("label")["dup_cluster"]
    m["dup_cluster"] = cl.loc[m["label"]].to_numpy()
    agg = {"disorder_class": ("disorder_class", "first"), "rg_ratio": ("rg_ratio", "median"),
           "chi2_raw": ("chi2_raw", "median")}
    collapsed = m.groupby("dup_cluster").agg(**agg).reset_index()
    out["near_duplicates"] = {
        "threshold_jaccard_5mer": robust.NEAR_DUPLICATE, "n_pairs": len(pairs),
        "n_identical_pairs": int(sum(p["identical"] for p in pairs)),
        "n_entries_in_clusters": int(m["dup_cluster"].duplicated(keep=False).sum()),
        "n_clusters_clean": int(collapsed.shape[0]),
        "pairs": pairs,
        "collapsed_curve_rg_ratio_by_class": _class_ratio_stats(collapsed, seed_offset=30),
        "collapsed_share_raw_fit": float((collapsed["chi2_raw"] <= 2).mean()),
        "collapsed_chi2_raw_median": float(collapsed["chi2_raw"].median())}
    p = C.RESULTS / "coord_rg.csv"
    if p.exists():
        cr = pd.read_csv(p)
        # the class statistics use the analysed sample only; worked examples are added for the figure
        cr = cr[cr["label"].isin(m["label"]) & (cr.get("reason", "") != "worked example")]
        rec = {"n": int(len(cr)), "by_class": {}}
        for i, cls in enumerate(CLASS_ORDER):
            g = cr[cr["disorder_class"] == cls]
            if len(g) < 3:
                continue
            rec["by_class"][cls] = {
                "n": int(len(g)),
                "coord_over_exp_median": float(g["ratio_coord_exp"].median()),
                "coord_over_exp_ci95": ci(g["ratio_coord_exp"], seed_offset=40 + i),
                "curve_over_exp_median": float(g["ratio_curve_exp"].median()),
                "curve_over_exp_ci95": ci(g["ratio_curve_exp"], seed_offset=50 + i),
                "curve_over_coord_median": float(g["ratio_curve_coord"].median()),
                "curve_over_coord_ci95": ci(g["ratio_curve_coord"], seed_offset=60 + i),
                "p_two_sided_wilcoxon_log_coord_over_exp":
                    float(stats.wilcoxon(np.log(g["ratio_coord_exp"])).pvalue) if len(g) > 5 else None}
        out["coordinate_rg"] = rec
    return out


def headline_robustness(d: pd.DataFrame) -> dict:
    """C1: bootstrap intervals for the raw-fit share overall and per class; C2: the error-free
    fit score and whether the class and model orderings hold under it; C3: kind counts against the
    phi threshold."""
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    fit = (m["chi2_raw"] <= 2).astype(float)
    out = {"share_raw_fit": {"value": float(fit.mean()), "ci95": ci(fit, robust.share, seed_offset=70),
                             "n": int(len(m))}}
    out["share_raw_fit_by_class"] = {
        cls: {"value": float((g["chi2_raw"] <= 2).mean()), "n": int(len(g)),
              "ci95": ci((g["chi2_raw"] <= 2).astype(float), robust.share, seed_offset=80 + i)}
        for i, (cls, g) in enumerate(m.groupby("disorder_class"))}
    out["chi2_raw_median_ci95"] = ci(m["chi2_raw"], seed_offset=90)
    if "nrmsd_raw" in m:
        by_cls = {cls: {"n": int(len(g)), "nrmsd_raw_median": float(g["nrmsd_raw"].median()),
                        "nrmsd_raw_ci95": ci(g["nrmsd_raw"], seed_offset=300 + i),
                        "nrmsd_best_median": float(g["nrmsd_best"].median()),
                        "chi2_raw_median": float(g["chi2_raw"].median())}
                  for i, (cls, g) in enumerate(m.groupby("disorder_class"))}
        order_chi2 = sorted(CLASS_ORDER, key=lambda c: by_cls[c]["chi2_raw_median"])
        order_nrmsd = sorted(CLASS_ORDER, key=lambda c: by_cls[c]["nrmsd_raw_median"])
        rho = stats.spearmanr(m["log_chi2_raw"], m["nrmsd_raw"])
        rank_chi2, rank_nrmsd = {}, {}
        base = m.set_index("label")
        for mo in sorted(d["model"].unique()):
            o = d[(d["model"] == mo) & d["clean"]].set_index("label")
            common = base.index.intersection(o.index)
            if len(common) < 10:
                continue
            rank_chi2[mo] = float(o.loc[common, "chi2_raw"].median())
            rank_nrmsd[mo] = float(o.loc[common, "nrmsd_raw"].median())
        out["nrmsd"] = {
            "definition": "root-mean-square deviation of ln I between the scaled model (chi2-fit scale and "
                          "background) and the data over points with I > 0, I >= 3 sigma and a positive "
                          "model, divided by the range of ln I of the data over those points",
            "by_class": by_cls, "class_order_chi2": order_chi2, "class_order_nrmsd": order_nrmsd,
            "class_order_holds": order_chi2 == order_nrmsd,
            "spearman_log_chi2_vs_nrmsd": float(rho.statistic),
            "model_median_chi2_raw": rank_chi2, "model_median_nrmsd_raw": rank_nrmsd,
            "model_order_chi2": sorted(rank_chi2, key=rank_chi2.get),
            "model_order_nrmsd": sorted(rank_nrmsd, key=rank_nrmsd.get),
            "model_order_holds": (sorted(rank_chi2, key=rank_chi2.get)
                                  == sorted(rank_nrmsd, key=rank_nrmsd.get)),
            "bioemu_first_under_nrmsd": min(rank_nrmsd, key=rank_nrmsd.get) == C.MODEL}
    out["kind_counts_by_phi_threshold"] = {
        str(t): m.apply(kind, axis=1, phi_threshold=t).value_counts().reindex(KINDS, fill_value=0).to_dict()
        for t in PHI_SENSITIVITY}
    return out


def main() -> None:
    ent, d = load()
    d["kind"] = d.apply(kind, axis=1)
    d["resolvability"] = d["kind"]
    clusters = sequence_clusters(ent)
    clusters.to_csv(C.RESULTS / "clusters.csv", index=False)
    bc = by_class(d)
    bc.to_csv(C.RESULTS / "by_class.csv", index=False)
    keep = ["model", "label", "disorder_class", "length", "chi2_raw", "chi2_best", "phi_at_chi2_1",
            "phi_at_chi2_2", "resolvability", "rg_exp", "rg_model", "rg_ratio", "clean"]
    keep += [c for c in ("nrmsd_raw", "nrmsd_best") if c in d]
    d[keep].to_csv(C.RESULTS / "resolvability.csv", index=False)
    cand = candidates(d)
    cand.to_csv(C.RESULTS / "candidates.csv", index=False)
    m = d[(d["model"] == C.MODEL)]
    mc = m[m["clean"]]
    summary = {
        "n_entries": int(len(ent)), "n_clean": int(clean_mask(ent).sum()),
        "n_flag_aggregation": int(ent["flag_aggregation"].sum()),
        "n_flag_negative": int(ent["flag_negative_I"].sum()),
        "n_flag_bad_errors": int(ent["flag_bad_errors"].sum()),
        "n_flag_oligomer": int(ent["flag_oligomer"].sum()),
        "oligomer_median_chi2_raw": float(m.loc[m["flag_oligomer"], "chi2_raw"].median()),
        "n_guinier_invalid": int((~ent["guinier_valid"]).sum()),
        "class_counts": ent["disorder_class"].value_counts().to_dict(),
        "bioemu": {
            "n_scored": int(len(m)), "n_conformers_median": float(m["n_conformers"].median()),
            "chi2_raw_median": float(mc["chi2_raw"].median()),
            "chi2_raw_q25": float(mc["chi2_raw"].quantile(0.25)),
            "chi2_raw_q75": float(mc["chi2_raw"].quantile(0.75)),
            "share_raw_fit_chi2_2": float((mc["chi2_raw"] <= 2).mean()),
            "share_raw_fit_chi2_1": float((mc["chi2_raw"] <= 1).mean()),
            "chi2_best_median": float(mc["chi2_best"].median()),
            "phi_at_chi2_1_median": float(mc["phi_at_chi2_1"].median()),
            "resolvability_counts": mc["kind"].value_counts().to_dict(),
            "rg_ratio_median_by_class": mc.groupby("disorder_class")["rg_ratio"].median().to_dict(),
        },
        "kinds": KINDS,
        "examples": choose_examples(d),
        "representative_path": choose_representative_path(d),
        "expectations": expectations(d),
        "rg_ratio_disordered_by_model": {
            m: float(g["rg_ratio"].median())
            for m, g in d[d["clean"] & (d["disorder_class"] == "disordered")].groupby("model")},
        "model_comparison": model_comparison(d),
        "by_class": json.loads(bc.to_json(orient="records")),
        "top_candidates": json.loads(cand.head(25).to_json(orient="records")),
        "thresholds": {"raw_fit_chi2": 2.0, "resolvability_phi": PHI_THRESHOLD, "chi2_target": C.CHI2_TARGET,
                       "aggregation_upturn": C.AGGREGATION_UPTURN,
                       "theta_grid": [float(x) for x in C.THETA_GRID]},
        "rg_robustness": rg_robustness(d, ent, clusters),
        "headline_robustness": headline_robustness(d),
    }
    json.dump(summary, open(C.RESULTS / "analysis.json", "w"), indent=1, default=float)
    b = summary["bioemu"]
    print(f"  {summary['n_clean']} clean entries; BioEmu raw chi2 median {b['chi2_raw_median']:.2f}; "
          f"raw fit (chi2<=2) in {b['share_raw_fit_chi2_2']:.0%}; "
          f"expectations met: {[k for k, v in summary['expectations'].items() if v['met']]}", flush=True)


if __name__ == "__main__":
    main()
