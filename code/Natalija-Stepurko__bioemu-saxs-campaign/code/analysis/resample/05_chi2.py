"""Raw χ² of the ensemble-average curve against the measured profile, by three paths and two conformer sets
(design §4 and §5).

    python 05_chi2.py

Paths: A, v1.2 through PULCHRA and Pepsi-SAXS (RESAMPLE_WORK/saxs/v12); B, the archive's conformers through
the same steps (RESAMPLE_WORK/saxs/archive); C, the archive's own computed curves (BSC_ROOT,
Pepsi-<label>.csv).
Sets: all conformers; the conformers that pass BioEmu's physical filter (results/resample/filter.csv; A uses
the v1.2 mask, B and C the archive mask). Scored with rslib.raw_fit, the raw fit of bsc.score.score_entry.

Checks, written with the results: path C over all conformers against the raw χ² of results/scores.csv; the
row order of the archive's curves against its conformers (Spearman correlation, per protein, of the Guinier
Rg of each row's curve with the Cα Rg of the conformer of the same index, for paths B and C).
Outputs: results/resample/chi2.csv, results/resample/curve_order.csv.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import rslib
from common import ARCHIVE_CURVES, OUT, RES, SAXS, labels, pin, record
from scipy import stats

from bsc import config as C
from bsc import saxs
from bsc.score import read_experiment

PATHS = {"A": ("v12", lambda lab: SAXS / "v12" / lab / "curves.csv"),
         "B": ("archive", lambda lab: SAXS / "archive" / lab / "curves.csv"),
         "C": ("archive", lambda lab: ARCHIVE_CURVES / f"Pepsi-{lab}.csv")}


def curve_rg(curves: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Guinier Rg of each conformer's curve, as bsc.score computes the per-conformer Rg."""
    return np.array([saxs.guinier(q, c, 0.01 * c + 1e-9, C.GUINIER_QRG_MAX)["rg"] if np.isfinite(c).all()
                     else np.nan for c in curves])


def main() -> None:
    t0 = time.time()
    used = pin()
    filt = pd.read_csv(OUT / "filter.csv")
    rgc = pd.read_csv(OUT / "rg_conformers.csv")
    rows, order = [], []
    for lab in labels():
        exp = read_experiment(lab)
        for path, (side, f) in PATHS.items():
            curves = pd.read_csv(f(lab), index_col=0).to_numpy(float)
            m = filt[(filt["label"] == lab) & (filt["side"] == side)].sort_values("conformer")
            if len(m) != len(curves):
                raise ValueError(f"{lab} path {path}: {len(curves)} curves, {len(m)} conformers")
            for subset, keep in (("all", np.ones(len(m), bool)), ("filtered", m["pass_all"].to_numpy())):
                r = rslib.raw_fit(curves[keep], exp) if keep.any() else {"n_conformers": 0, "n_dropped": 0,
                                                                         "chi2": np.nan, "nrmsd": np.nan}
                rows.append({"label": lab, "path": path, "side": side, "subset": subset, **r})
            if side == "archive":
                g = rgc[(rgc["label"] == lab) & (rgc["side"] == side)].sort_values("conformer")
                g = g["rg_ca"].to_numpy()
                crg = curve_rg(curves, exp["q"].to_numpy())
                ok = np.isfinite(crg)
                order.append({"label": lab, "path": path, "n": int(ok.sum()),
                              "spearman_curve_rg_vs_ca_rg": float(stats.spearmanr(crg[ok], g[ok])[0])})
        print(lab, flush=True)
    out = pd.DataFrame(rows)
    sc = pd.read_csv(RES / "scores.csv")
    sc = sc[sc["model"] == "bioemu"].set_index("label")["chi2_raw"]
    c_all = out[(out["path"] == "C") & (out["subset"] == "all")].set_index("label")
    rel = (c_all["chi2"] / sc.loc[c_all.index] - 1).abs()
    out.to_csv(OUT / "chi2.csv", index=False)
    pd.DataFrame(order).to_csv(OUT / "curve_order.csv", index=False)
    record("05_chi2", {"cpus": used, "wall_seconds": round(time.time() - t0, 1),
                       "path_C_vs_scores_csv_max_relative_difference": float(rel.max())})
    if not rel.max() < 1e-6:
        sys.exit(f"path C does not reproduce results/scores.csv (relative difference up to {rel.max():.1e})")


if __name__ == "__main__":
    main()
