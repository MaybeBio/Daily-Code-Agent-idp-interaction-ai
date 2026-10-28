"""Analysis and export of the v1.2 re-sample (design §3–§6).

    python 06_analyse.py

Reads results/resample/{inventory,conversion,filter,rg_conformers,chi2,curve_order}.csv, the tracked
results/coord_rg.csv and results/resample/versions.json. Writes results/resample/per_protein.csv (one row
per protein), results/resample/summary.json (every number reported) and results/figures/fig_resample_rg.png.
"""

from __future__ import annotations

import json
import time

import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import rslib
from common import FIG, MIN_FILTERED, N_BOOT, OUT, RES, SEED, VERSIONS, record
from scipy import stats

from bsc import report

CLASSES = ("disordered", "folded")
LENGTH_CUT = 100      # residues; the length step the study reports for folded proteins (E2), used to describe
                      # the disordered ratios by length (added after the first results, design §9)


def ensemble_table(rgc: pd.DataFrame, filt: pd.DataFrame) -> pd.DataFrame:
    """Per protein and side: Rg_ens over all conformers and over those passing the filter."""
    m = rgc.merge(filt[["label", "side", "conformer", "pass_all"]], on=["label", "side", "conformer"],
                  how="left", validate="one_to_one")
    if m["pass_all"].isna().any():
        raise ValueError("conformers without a filter result")
    rows = []
    for (lab, side), g in m.groupby(["label", "side"], sort=False):
        k = g["pass_all"].astype(bool).to_numpy()
        rows.append({"label": lab, "side": side, "n": len(g), "n_pass": int(k.sum()),
                     "rg_ens_all": rslib.ensemble_rg(g["rg_ca"].to_numpy()),
                     "rg_ens_filtered": rslib.ensemble_rg(g["rg_ca"].to_numpy()[k])})
    return pd.DataFrame(rows)


def build() -> pd.DataFrame:
    inv = pd.read_csv(OUT / "inventory.csv")
    conv = pd.read_csv(OUT / "conversion.csv")
    filt = pd.read_csv(OUT / "filter.csv")
    rgc = pd.read_csv(OUT / "rg_conformers.csv")
    chi = pd.read_csv(OUT / "chi2.csv")
    ens = ensemble_table(rgc, filt)
    d = inv[["label", "length", "disorder_class", "rg_exp", "guinier_valid", "colab_seconds",
             "has_samples_xtc", "has_unphysical_xtc", "frames_samples_xtc", "frames_unphysical_xtc"]].copy()
    d = d.merge(conv[["label", "n_recovered"]], on="label")
    for side, tag in (("v12", "v12"), ("archive", "arch")):
        e = ens[ens["side"] == side].set_index("label")
        f = filt[filt["side"] == side]
        d[f"n_{tag}"] = e.loc[d["label"], "n"].to_numpy()
        d[f"n_pass_{tag}"] = e.loc[d["label"], "n_pass"].to_numpy()
        d[f"rejected_share_{tag}"] = 1 - d[f"n_pass_{tag}"] / d[f"n_{tag}"]
        for crit in ("ca_seq", "cn_seq", "clash"):
            fail = (~f.groupby("label")[f"pass_{crit}"].all()).rename(None)
            share = 1 - f.groupby("label")[f"pass_{crit}"].mean()
            d[f"fail_{crit}_share_{tag}"] = share.loc[d["label"]].to_numpy()
            d[f"any_fail_{crit}_{tag}"] = fail.loc[d["label"]].to_numpy()
        for sub in ("all", "filtered"):
            d[f"rg_ca_{tag}_{sub}"] = e.loc[d["label"], f"rg_ens_{sub}"].to_numpy()
            d[f"ratio_{tag}_{sub}"] = d[f"rg_ca_{tag}_{sub}"] / d["rg_exp"]
    for path in ("A", "B", "C"):
        for sub in ("all", "filtered"):
            c = chi[(chi["path"] == path) & (chi["subset"] == sub)].set_index("label")
            for col in ("chi2", "nrmsd", "n_conformers", "n_dropped"):
                d[f"{col}_{path}_{sub}"] = c.loc[d["label"], col].to_numpy()
    d["filtered_comparable"] = (d["n_pass_v12"] >= MIN_FILTERED) & (d["n_pass_arch"] >= MIN_FILTERED)
    cr = pd.read_csv(RES / "coord_rg.csv").set_index("label")["rg_ca_rms"]
    d["rg_ca_arch_coord_rg_csv"] = cr.reindex(d["label"]).to_numpy()
    return d


def by_class_shares(d: pd.DataFrame, tag: str) -> dict:
    out = {}
    for cls, g in d.groupby("disorder_class"):
        s = g[f"rejected_share_{tag}"]
        out[cls] = {"n_proteins": len(g), "median": float(s.median()), "mean": float(s.mean()),
                    "min": float(s.min()), "max": float(s.max()),
                    "n_with_any_rejected": int((s > 0).sum()), "n_all_rejected": int((s == 1).sum()),
                    "mean_fail_share_by_criterion": {c: float(g[f"fail_{c}_share_{tag}"].mean())
                                                     for c in ("ca_seq", "cn_seq", "clash")},
                    "n_proteins_with_any_failure_by_criterion": {c: int(g[f"any_fail_{c}_{tag}"].sum())
                                                                 for c in ("ca_seq", "cn_seq", "clash")}}
    return out


def size(d: pd.DataFrame, sub: str) -> dict:
    out = {}
    for cls in CLASSES:
        g = d[d["disorder_class"] == cls]
        if sub == "filtered":
            g = g[g["filtered_comparable"]]
        out[cls] = rslib.paired(g[f"ratio_v12_{sub}"], g[f"ratio_arch_{sub}"], SEED, N_BOOT)
    for _, r in d[~d["disorder_class"].isin(CLASSES)].iterrows():
        out[f"{r['disorder_class']} ({r['label']})"] = {"ratio_v12": r[f"ratio_v12_{sub}"],
                                                         "ratio_archive": r[f"ratio_arch_{sub}"]}
    return out


def length_split(d: pd.DataFrame) -> dict:
    """Disordered proteins by chain length: the Cα-Rg ratios of both sides and v1.2's rejected share, for
    chains of at most LENGTH_CUT residues and longer ones, and the rank correlation of each with length."""
    g = d[d["disorder_class"] == "disordered"]
    out = {"cut_residues": LENGTH_CUT}
    for name, h in (("short", g[g["length"] <= LENGTH_CUT]), ("long", g[g["length"] > LENGTH_CUT])):
        out[name] = {"n": len(h), "length_min": int(h["length"].min()), "length_max": int(h["length"].max()),
                     "ratio_v12": rslib.median_ci(h["ratio_v12_all"], SEED, N_BOOT),
                     "ratio_archive": rslib.median_ci(h["ratio_arch_all"], SEED, N_BOOT),
                     "n_v12_below_1": int((h["ratio_v12_all"] < 1).sum()),
                     "n_v12_below_archive": int((h["ratio_v12_all"] < h["ratio_arch_all"]).sum()),
                     "rejected_share_v12_median": float(h["rejected_share_v12"].median()),
                     "n_all_rejected_v12": int((h["n_pass_v12"] == 0).sum())}
    for key, col in (("ratio_v12", "ratio_v12_all"), ("ratio_archive", "ratio_arch_all"),
                     ("rejected_share_v12", "rejected_share_v12")):
        r = stats.spearmanr(g["length"], g[col])
        out[f"spearman_length_vs_{key}"] = {"rho": float(r.statistic), "p": float(r.pvalue), "n": len(g)}
    return out


def fit(d: pd.DataFrame, sub: str) -> dict:
    g0 = d[d["filtered_comparable"]] if sub == "filtered" else d
    out = {}
    for name, g in [("all proteins", g0)] + [(c, g0[g0["disorder_class"] == c]) for c in CLASSES]:
        out[name] = {"A_vs_B": rslib.paired(g[f"chi2_A_{sub}"], g[f"chi2_B_{sub}"], SEED, N_BOOT),
                     "A_vs_C": rslib.paired(g[f"chi2_A_{sub}"], g[f"chi2_C_{sub}"], SEED, N_BOOT),
                     "B_vs_C": rslib.paired(g[f"chi2_B_{sub}"], g[f"chi2_C_{sub}"], SEED, N_BOOT),
                     "nrmsd_A_vs_B": rslib.paired(g[f"nrmsd_A_{sub}"], g[f"nrmsd_B_{sub}"], SEED, N_BOOT)}
    return out


def figure(d: pd.DataFrame) -> None:
    report.style()
    fig, ax = report.figure(560, 470)
    v = d[["ratio_arch_all", "ratio_v12_all"]].to_numpy()
    lo, hi = min(0.55, 0.92 * v.min()), max(2.4, 1.08 * v.max())     # every protein inside the axes
    ax.plot([lo, hi], [lo, hi], color=report.MUTED, lw=0.9, zorder=1, label="same size on both sides")
    ax.axhline(1, color=report.RULE, lw=0.8, ls=(0, (3, 3)), zorder=0)
    ax.axvline(1, color=report.RULE, lw=0.8, ls=(0, (3, 3)), zorder=0)
    for cls in report.CLASS_ORDER:
        g = d[d["disorder_class"] == cls]
        if g.empty:
            continue
        full = g["n_pass_v12"] == 0
        ax.scatter(g.loc[~full, "ratio_arch_all"], g.loc[~full, "ratio_v12_all"], s=34,
                   color=report.CLS_COL[cls], alpha=0.8, edgecolor="white", linewidth=0.6, zorder=3,
                   label=f"{cls} (n = {len(g)})")
        if full.any():
            ax.scatter(g.loc[full, "ratio_arch_all"], g.loc[full, "ratio_v12_all"], s=34, facecolor="white",
                       edgecolor=report.CLS_COL[cls], linewidth=1.4, zorder=4)
    ax.scatter([], [], s=34, facecolor="white", edgecolor=report.MUTED, linewidth=1.4,
               label="every v1.2 sample rejected\nby BioEmu's filter")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = [t for t in (0.5, 0.6, 0.8, 1.0, 1.25, 1.6, 2.0, 2.5) if lo <= t <= hi]
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(mticker.FixedLocator(ticks))
        axis.set_major_formatter(mticker.FixedFormatter([f"{t:g}" for t in ticks]))
        axis.set_minor_locator(mticker.NullLocator())
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("archive ensemble: Cα Rg / measured Guinier Rg")
    ax.set_ylabel("bioemu-v1.2: Cα Rg / measured Guinier Rg")
    ax.legend(loc="upper left", frameon=False, handletextpad=0.3)
    fig.tight_layout()
    report.OUT = FIG
    report.save(fig, "fig_resample_rg")


def main() -> None:
    t0 = time.time()
    d = build()
    d.to_csv(OUT / "per_protein.csv", index=False)
    conv = pd.read_csv(OUT / "conversion.csv")
    order = pd.read_csv(OUT / "curve_order.csv")
    ver = json.loads(VERSIONS.read_text())
    phys = d[d["has_samples_xtc"]]
    unph = d[d["has_unphysical_xtc"]]
    full = d[d["n_pass_v12"] == 0].sort_values("length")
    rho = order.groupby("path")["spearman_curve_rg_vs_ca_rg"]
    S = {
        "proteins": {"n": len(d), "by_class": d["disorder_class"].value_counts().to_dict()},
        "files": {
            "n_with_samples_xtc": int(d["has_samples_xtc"].sum()),
            "n_with_only_unphysical_xtc": int((d["has_unphysical_xtc"] & ~d["has_samples_xtc"]).sum()),
            "n_with_both": int((d["has_unphysical_xtc"] & d["has_samples_xtc"]).sum()),
            "samples_xtc_frames_equal_filter_passes": bool(
                (phys["frames_samples_xtc"] == phys["n_pass_v12"]).all()),
            "unphysical_xtc_frames_equal_all_samples": bool(
                (unph["frames_unphysical_xtc"] == unph["n_recovered"]).all()),
            "unphysical_xtc_frames": sorted(unph["frames_unphysical_xtc"].astype(int).unique().tolist()),
            "reproduction_files_match": bool(conv["files_match"].all()),
            "reproduction_max_abs_coordinate_difference_A": float(np.nanmax(
                conv.filter(like="max_abs_diff_A_").to_numpy(float))),
            "max_rmsd_A_colab_frame_vs_recovered_sample": float(np.nanmax(
                conv.filter(like="max_rmsd_A_colab_vs_recovered_").to_numpy(float))),
            "max_mdtraj_float32_rmsd_A_same_comparison": float(np.nanmax(
                conv.filter(like="max_mdtraj_rmsd_A_colab_vs_recovered_").to_numpy(float))),
            "n_proteins_with_rmsd_check": int(conv.filter(like="max_rmsd_A_colab_vs_recovered_").notna()
                                              .any(axis=1).sum()),
        },
        "recovered": {"total": int(d["n_recovered"].sum()), "per_protein_min": int(d["n_recovered"].min()),
                      "per_protein_max": int(d["n_recovered"].max()),
                      "archive_conformers_min": int(d["n_arch"].min()),
                      "archive_conformers_max": int(d["n_arch"].max()),
                      "kept_by_filter_on_colab_total": int(phys["frames_samples_xtc"].sum())},
        "rejected_v12": by_class_shares(d, "v12"),
        "rejected_archive_same_criteria": by_class_shares(d, "arch"),
        "all_rejected_v12": [{"label": r["label"], "length": int(r["length"]), "class": r["disorder_class"],
                              **{f"fail_{c}_share": float(r[f"fail_{c}_share_v12"])
                                 for c in ("ca_seq", "cn_seq", "clash")},
                              "archive_rejected_share": float(r["rejected_share_arch"])}
                             for _, r in full.iterrows()],
        "all_rejected_length_vs_rest": {"median_length_all_rejected": float(full["length"].median()),
                                        "median_length_rest": float(
                                            d.loc[d["n_pass_v12"] > 0, "length"].median())},
        "size_all": size(d, "all"),
        "size_filtered": size(d, "filtered"),
        # all samples, on the proteins of the filtered comparison: separates the filter from the protein set
        "size_all_filtered_set": size(d[d["filtered_comparable"]], "all"),
        "length_split_disordered": length_split(d),
        "n_filtered_comparable": {"total": int(d["filtered_comparable"].sum()),
                                  "by_class": d[d["filtered_comparable"]]["disorder_class"]
                                  .value_counts().to_dict()},
        "fit_all": fit(d, "all"),
        "fit_filtered": fit(d, "filtered"),
        "fit_all_filtered_set": fit(d[d["filtered_comparable"]], "all"),
        "dropped_curves": {p: int(d[f"n_dropped_{p}_all"].sum()) for p in ("A", "B", "C")},
        "checks": {
            "archive_rg_vs_coord_rg_csv_n": int(d["rg_ca_arch_coord_rg_csv"].notna().sum()),
            "archive_rg_vs_coord_rg_csv_max_abs_diff_A": float(
                (d["rg_ca_arch_all"] - d["rg_ca_arch_coord_rg_csv"]).abs().max()),
            "path_C_vs_scores_csv_max_relative_difference":
                ver["05_chi2"]["path_C_vs_scores_csv_max_relative_difference"],
            "curve_order_spearman_median": rho.median().to_dict(),
            "curve_order_spearman_min": rho.min().to_dict(),
            "n_guinier_valid": int(d["guinier_valid"].sum()),
        },
        "timings_seconds": {
            "colab_sampling_total": float(d["colab_seconds"].sum()),
            "colab_sampling_median_per_protein": float(d["colab_seconds"].median()),
            "colab_sampling_max_per_protein": float(d["colab_seconds"].max()),
            **{k: v.get("wall_seconds", v.get("wall_seconds_this_run")) for k, v in ver.items()
               if k != "06_analyse"},
            "04_saxs_seconds_per_task_sum": ver.get("04_saxs", {}).get("seconds_per_task_sum"),
        },
    }
    e12 = S["size_all"]["disordered"]["a"]
    S["E12"] = {"statement": "v1.2's median Cα-Rg ratio over the disordered proteins, all samples, "
                             "is above 1.1",
                "median": e12["median"], "ci95": e12["ci95"], "n": e12["n"], "met": bool(e12["median"] > 1.1),
                "lower_bound_above_1_1": bool(e12["ci95"][0] > 1.1)}
    e13 = S["fit_all"]["all proteins"]["A_vs_B"]
    S["E13"] = {"statement": "v1.2's median raw χ² over all listed proteins (path A) is lower than the "
                             "archive's through the same path (B)",
                "median_A": e13["a"]["median"], "median_B": e13["b"]["median"], "n": e13["n"],
                "median_A_minus_median_B": e13["median_a_minus_median_b"],
                "ci95": e13["median_a_minus_median_b_ci95"],
                "wilcoxon_p": e13.get("wilcoxon_two_sided", {}).get("p"),
                "met": bool(e13["a"]["median"] < e13["b"]["median"])}
    (OUT / "summary.json").write_text(json.dumps(S, indent=1, default=float) + "\n")
    figure(d)
    record("06_analyse", {"seed": SEED, "n_boot": N_BOOT, "min_filtered": MIN_FILTERED,
                          "wall_seconds": round(time.time() - t0, 1)})
    print(json.dumps({"E12": S["E12"], "E13": S["E13"]}, indent=1))


if __name__ == "__main__":
    main()
