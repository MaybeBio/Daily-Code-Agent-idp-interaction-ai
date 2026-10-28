"""Stage `score`: per entry, the raw-ensemble fit, the Guinier comparison, the reweighting path and
the protein class, for the model under study and each comparator present in the archive.

Output: results/scores.csv (one row per entry × model), results/paths/<model>.json (the full
reweighting path per entry), results/entries.csv (metadata and data-quality flags per entry).
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from bsc import config as C
from bsc import saxs


def read_experiment(label: str) -> pd.DataFrame:
    """q, I(q), sigma of a cleaned SASBDB curve (PeptoneDB format: a comment line, then three columns)."""
    df = pd.read_csv(C.SAXS_DIR / f"{label}-bift.dat", sep=r"\s+", comment="#", header=None,
                     names=["q", "I", "sigma"])
    return df


def read_prediction(model: str, label: str) -> np.ndarray | None:
    """Back-calculated curves, shape (n_conformers, n_q); columns are q values matching the data."""
    p = C.PRED_DIR / model / f"{C.PREDICTOR}-{label}.csv"
    if not p.exists():
        return None
    return pd.read_csv(p, index_col=0).to_numpy(float)


AA_MASS = {"A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14, "E": 129.12, "Q": 128.13,
           "G": 57.05, "H": 137.14, "I": 113.16, "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18,
           "P": 97.12, "S": 87.08, "T": 101.10, "W": 186.21, "Y": 163.18, "V": 99.13}


def sequence_mass_kda(seq: str) -> float:
    return (sum(AA_MASS.get(a, 110.0) for a in seq) + 18.02) / 1000


def disorder_class(mean_gscore: float) -> str:
    for name, lo, hi in C.DISORDER_CLASSES:
        if lo <= mean_gscore < hi:
            return name
    return "unknown"


def entry_metadata(table: pd.DataFrame) -> pd.DataFrame:
    summ = json.load(open(C.SASBDB_SUMMARY)) if C.SASBDB_SUMMARY.exists() else {}
    rows = []
    for _, r in table.iterrows():
        sm = summ.get(r["label"], {})
        mw_exp = float(sm["experimental_mw"]) if sm.get("experimental_mw") else np.nan
        mw_seq = sequence_mass_kda(r["sequence"])
        exp = read_experiment(r["label"])
        usable = exp[exp["sigma"] > 0]
        g = (saxs.guinier(usable["q"].to_numpy(), usable["I"].to_numpy(), usable["sigma"].to_numpy(),
                          C.GUINIER_QRG_MAX) if (usable["I"] > 0).sum() >= 12
             else {"rg": np.nan, "rg_err": np.nan, "n_points": 0, "upturn": np.nan, "valid": False})
        rows.append({"label": r["label"], "length": int(r["length"]), "ph": r.get("pH", np.nan),
                     "disorder_mean": float(r["mean_gscore_adopt2"]),
                     "disorder_class": disorder_class(float(r["mean_gscore_adopt2"])),
                     "length_bin": str(pd.cut([r["length"]], C.LENGTH_BINS).astype(str)[0]),
                     "n_q": len(exp), "q_min": float(exp["q"].min()), "q_max": float(exp["q"].max()),
                     "rg_exp": g["rg"], "rg_exp_err": g["rg_err"], "guinier_points": g["n_points"],
                     "guinier_valid": g["valid"], "upturn": g["upturn"],
                     "flag_aggregation": bool(np.isfinite(g["upturn"])
                                              and g["upturn"] > C.AGGREGATION_UPTURN),
                     "flag_negative_I": bool(((usable["I"] / usable["sigma"]) < -C.NEGATIVE_SIGMA).any()),
                     "flag_bad_errors": bool((exp["sigma"] <= 0).any()),
                     "mw_seq_kda": mw_seq, "mw_exp_kda": mw_exp, "mw_ratio": mw_exp / mw_seq,
                     "flag_oligomer": bool(np.isfinite(mw_exp) and mw_exp / mw_seq > C.OLIGOMER_MW_RATIO),
                     "sasbdb_title": sm.get("title"),
                     "sequence": r["sequence"]})
    return pd.DataFrame(rows)


def score_entry(model: str, label: str, exp: pd.DataFrame) -> tuple[dict, list[dict]] | None:
    curves = read_prediction(model, label)
    if curves is None:
        return None
    if (exp["sigma"] <= 0).all():
        return None
    q, I, s = exp["q"].to_numpy(), exp["I"].to_numpy(), exp["sigma"].to_numpy()
    s = np.where(s > 0, s, np.nan)
    keep = np.isfinite(s)
    if curves.shape[1] != len(q):
        raise ValueError(f"{model} {label}: {curves.shape[1]} q points in the prediction, "
                         f"{len(q)} in the data")
    ok = np.isfinite(curves).all(axis=1)
    curves = curves[ok][:, keep]
    q, I, s = q[keep], I[keep], s[keep]
    raw = saxs.ensemble_curve(curves)
    rec = {"model": model, "label": label, "n_conformers": int(ok.sum()), "n_dropped": int((~ok).sum()),
           "chi2_raw": saxs.chi2(raw, I, s),
           "chi2_raw_noback": saxs.chi2(raw, I, s, fit_background=False)}
    # Guinier Rg of the raw ensemble curve, on the experimental q grid
    g = saxs.guinier(q, raw, 0.01 * raw + 1e-9, C.GUINIER_QRG_MAX)
    rec["rg_model"] = g["rg"]
    # per-conformer Rg from each conformer's own curve
    rgs = [saxs.guinier(q, c, 0.01 * c + 1e-9, C.GUINIER_QRG_MAX)["rg"] for c in curves]
    rec["rg_conformer_median"] = float(np.nanmedian(rgs))
    rec["rg_conformer_iqr"] = float(np.nanpercentile(rgs, 75) - np.nanpercentile(rgs, 25))
    path, weights = saxs.reweighting_curve(curves, I, s, C.THETA_GRID, return_weights=True)
    i_best = int(np.argmin([p["chi2"] for p in path]))
    best = path[i_best]
    # a second fit score that does not weight by the reported errors (DESIGN §7, C2)
    rec["nrmsd_raw"] = saxs.nrmsd_log(raw, I, s)
    rec["nrmsd_best"] = saxs.nrmsd_log(saxs.ensemble_curve(curves, weights[i_best]), I, s)
    rec.update({"chi2_best": best["chi2"], "phi_at_best": best["phi"],
                "phi_at_chi2_1": saxs.phi_at_chi2(path, C.CHI2_TARGET),
                "phi_at_chi2_2": saxs.phi_at_chi2(path, 2.0),
                "chi2_at_phi_0.5": saxs.chi2_at_phi(path, 0.5),
                "chi2_at_phi_0.1": saxs.chi2_at_phi(path, 0.1)})
    # the benchmark's own operating point: the weights with effective sample size ESS_TARGET
    phi_bench = C.ESS_TARGET_BENCHMARK / len(curves)
    rec["chi2_at_phi_bench"] = saxs.chi2_at_phi(path, phi_bench)
    return rec, [{"label": label, **p} for p in path]


def main(models: list[str] | None = None) -> None:
    C.RESULTS.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(C.SAXS_TABLE)
    meta = entry_metadata(table)
    meta.to_csv(C.RESULTS / "entries.csv", index=False)
    print(f"  {len(meta)} entries; {meta['flag_aggregation'].sum()} flagged for aggregation, "
          f"{(~meta['guinier_valid']).sum()} without a valid Guinier region", flush=True)
    exps = {lab: read_experiment(lab) for lab in meta["label"]}
    models = models or [m for m in [C.MODEL, *C.COMPARATORS] if (C.PRED_DIR / m).exists()]
    rows = []
    (C.RESULTS / "paths").mkdir(exist_ok=True)
    for m in models:
        out = Parallel(n_jobs=C.N_JOBS)(delayed(score_entry)(m, lab, exps[lab]) for lab in meta["label"])
        out = [o for o in out if o is not None]
        rows += [r for r, _ in out]
        json.dump([p for _, path in out for p in path], open(C.RESULTS / "paths" / f"{m}.json", "w"))
        df = pd.DataFrame([r for r, _ in out])
        print(f"  {m}: {len(df)} entries; median raw chi2 {df['chi2_raw'].median():.2f}, "
              f"median best chi2 {df['chi2_best'].median():.2f}", flush=True)
    pd.DataFrame(rows).to_csv(C.RESULTS / "scores.csv", index=False)


if __name__ == "__main__":
    main()
